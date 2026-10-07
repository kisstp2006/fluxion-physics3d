// SPDX-License-Identifier: BSD-2-Clause

//! A tree of boxes, for finding what overlaps a box or a ray without
//! looking at everything: the world keeps one for what never moves and one
//! for what does, and a mesh keeps one over its triangles.
//!
//! ```zig
//! var tree: Tree = .empty;
//! defer tree.deinit(gpa);
//! const proxy = try tree.insert(gpa, wall_box, wall_index);
//! tree.query(player_box, &hits, Hits.add);
//! tree.remove(proxy);
//! ```
//!
//! **A bounding volume hierarchy**: every leaf is one thing's box, and every
//! branch holds the box round both its children. Asking what overlaps a box
//! walks down only the branches whose boxes do; a ray walks down only the
//! branches it crosses.
//!
//! **Built one leaf at a time, and kept balanced** - Erin Catto's dynamic
//! tree. A new leaf goes down towards the sibling whose box would grow least
//! by taking it in, measured by surface area: the chance that a random ray
//! or box meets a region is in proportion to its surface. Walking back up,
//! a branch whose two children differ in height by more than one is
//! rotated, as an AVL tree rotates, so no order of insertion grows it into a
//! list.
//!
//! Leaves stay where they are when others come and go, so a leaf's index -
//! its *proxy* - names it until it is removed.

const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Aabb = geometry.Aabb;

pub const Tree = @This();

/// No node: an empty tree's root, a leaf's missing children, the end of the
/// free list.
pub const null_node = std.math.maxInt(u32);

nodes: std.ArrayList(Node) = .empty,
root: u32 = null_node,
/// The first of the nodes not in use, linked through `parent`.
free: u32 = null_node,
leaf_count: u32 = 0,

pub const empty: Tree = .{};

const Node = struct {
    box: Aabb,
    /// The parent, or while the node is free, the next free node.
    parent: u32,
    child1: u32 = null_node,
    child2: u32 = null_node,
    /// Zero for a leaf, one more than its taller child for a branch, minus
    /// one while free.
    height: i32 = 0,
    /// What a leaf stands for: the caller's index.
    value: u32 = 0,

    fn isLeaf(self: *const Node) bool {
        return self.child1 == null_node;
    }
};

pub fn deinit(self: *Tree, gpa: Allocator) void {
    self.nodes.deinit(gpa);
    self.* = undefined;
}

/// Everything out, keeping the memory.
pub fn clear(self: *Tree) void {
    self.nodes.clearRetainingCapacity();
    self.root = null_node;
    self.free = null_node;
    self.leaf_count = 0;
}

// -------------------------------------------------------------------------
// Changing it
// -------------------------------------------------------------------------

/// Put a box in, standing for `value`, and hand back its proxy.
pub fn insert(self: *Tree, gpa: Allocator, box: Aabb, value: u32) Allocator.Error!u32 {
    const leaf = try self.allocate(gpa);
    self.nodes.items[leaf] = .{ .box = box, .parent = null_node, .value = value };
    try self.insertLeaf(gpa, leaf);
    self.leaf_count += 1;
    return leaf;
}

/// Take a leaf out. Its proxy means nothing afterwards.
pub fn remove(self: *Tree, proxy: u32) void {
    self.removeLeaf(proxy);
    self.release(proxy);
    self.leaf_count -= 1;
}

/// Give a leaf a new box: out and back in, keeping its proxy.
pub fn move(self: *Tree, gpa: Allocator, proxy: u32, box: Aabb) Allocator.Error!void {
    self.removeLeaf(proxy);
    self.nodes.items[proxy].box = box;
    try self.insertLeaf(gpa, proxy);
}

pub fn boxOf(self: *const Tree, proxy: u32) Aabb {
    return self.nodes.items[proxy].box;
}

pub fn valueOf(self: *const Tree, proxy: u32) u32 {
    return self.nodes.items[proxy].value;
}

/// How deep the tree is: zero for one leaf.
pub fn height(self: *const Tree) i32 {
    if (self.root == null_node) return 0;
    return self.nodes.items[self.root].height;
}

fn allocate(self: *Tree, gpa: Allocator) Allocator.Error!u32 {
    if (self.free != null_node) {
        const index = self.free;
        self.free = self.nodes.items[index].parent;
        return index;
    }
    const index: u32 = @intCast(self.nodes.items.len);
    try self.nodes.append(gpa, .{ .box = undefined, .parent = null_node });
    return index;
}

fn release(self: *Tree, index: u32) void {
    self.nodes.items[index] = .{ .box = undefined, .parent = self.free, .height = -1 };
    self.free = index;
}

fn cost(box: Aabb) f32 {
    return box.surfaceArea();
}

fn insertLeaf(self: *Tree, gpa: Allocator, leaf: u32) Allocator.Error!void {
    if (self.root == null_node) {
        self.root = leaf;
        self.nodes.items[leaf].parent = null_node;
        return;
    }

    // Down to the best sibling. At each branch: a new parent for this branch
    // and the leaf here, or on down into one child or the other - which
    // enlarges this branch by the leaf either way.
    const box = self.nodes.items[leaf].box;
    var index = self.root;
    while (!self.nodes.items[index].isLeaf()) {
        const node = &self.nodes.items[index];
        const combined = cost(node.box.join(box));
        const here = 2 * combined;
        const inherited = 2 * (combined - cost(node.box));
        const cost1 = self.descentCost(node.child1, box) + inherited;
        const cost2 = self.descentCost(node.child2, box) + inherited;
        if (here < cost1 and here < cost2) break;
        index = if (cost1 < cost2) node.child1 else node.child2;
    }
    const sibling = index;

    // `allocate` may move the nodes, so no pointer into them lives across it.
    const old_parent = self.nodes.items[sibling].parent;
    const branch = try self.allocate(gpa);
    const nodes = self.nodes.items;
    nodes[branch] = .{
        .box = box.join(nodes[sibling].box),
        .parent = old_parent,
        .child1 = sibling,
        .child2 = leaf,
        .height = nodes[sibling].height + 1,
    };
    nodes[sibling].parent = branch;
    nodes[leaf].parent = branch;
    if (old_parent == null_node) {
        self.root = branch;
    } else if (nodes[old_parent].child1 == sibling) {
        nodes[old_parent].child1 = branch;
    } else {
        nodes[old_parent].child2 = branch;
    }

    self.refit(nodes[leaf].parent);
}

fn descentCost(self: *const Tree, child: u32, box: Aabb) f32 {
    const node = &self.nodes.items[child];
    const grown = cost(node.box.join(box));
    return if (node.isLeaf()) grown else grown - cost(node.box);
}

fn removeLeaf(self: *Tree, leaf: u32) void {
    if (leaf == self.root) {
        self.root = null_node;
        return;
    }
    const nodes = self.nodes.items;
    const parent = nodes[leaf].parent;
    const grandparent = nodes[parent].parent;
    const sibling = if (nodes[parent].child1 == leaf) nodes[parent].child2 else nodes[parent].child1;

    if (grandparent == null_node) {
        self.root = sibling;
        nodes[sibling].parent = null_node;
        self.release(parent);
        return;
    }
    if (nodes[grandparent].child1 == parent) nodes[grandparent].child1 = sibling else nodes[grandparent].child2 = sibling;
    nodes[sibling].parent = grandparent;
    self.release(parent);
    self.refit(grandparent);
}

fn refit(self: *Tree, start: u32) void {
    var index = start;
    while (index != null_node) {
        index = self.balance(index);
        const nodes = self.nodes.items;
        const a = nodes[index].child1;
        const b = nodes[index].child2;
        nodes[index].height = 1 + @max(nodes[a].height, nodes[b].height);
        nodes[index].box = nodes[a].box.join(nodes[b].box);
        index = nodes[index].parent;
    }
}

fn balance(self: *Tree, a: u32) u32 {
    const nodes = self.nodes.items;
    if (nodes[a].isLeaf() or nodes[a].height < 2) return a;
    const b = nodes[a].child1;
    const c = nodes[a].child2;
    const lean = nodes[c].height - nodes[b].height;
    if (lean > 1) return self.rotate(a, c, .right);
    if (lean < -1) return self.rotate(a, b, .left);
    return a;
}

const Side = enum { left, right };

/// Lift `up`, a child of `a`, into `a`'s place; its taller child stays
/// with it and the shorter goes down to `a`.
fn rotate(self: *Tree, a: u32, up: u32, side: Side) u32 {
    const nodes = self.nodes.items;
    const other = if (side == .right) nodes[a].child1 else nodes[a].child2;
    const f = nodes[up].child1;
    const g = nodes[up].child2;

    nodes[up].child1 = a;
    nodes[up].parent = nodes[a].parent;
    nodes[a].parent = up;
    if (nodes[up].parent == null_node) {
        self.root = up;
    } else if (nodes[nodes[up].parent].child1 == a) {
        nodes[nodes[up].parent].child1 = up;
    } else {
        nodes[nodes[up].parent].child2 = up;
    }

    const keep, const give = if (nodes[f].height > nodes[g].height) .{ f, g } else .{ g, f };
    nodes[up].child2 = keep;
    if (side == .right) nodes[a].child2 = give else nodes[a].child1 = give;
    nodes[give].parent = a;

    nodes[a].box = nodes[other].box.join(nodes[give].box);
    nodes[a].height = 1 + @max(nodes[other].height, nodes[give].height);
    nodes[up].box = nodes[a].box.join(nodes[keep].box);
    nodes[up].height = 1 + @max(nodes[a].height, nodes[keep].height);
    return up;
}

// -------------------------------------------------------------------------
// Asking it
// -------------------------------------------------------------------------

/// Deep enough for any tree balance allows.
const stack_depth = 128;

/// Call `visit(context, value)` for every leaf whose box overlaps `box`,
/// until it returns false.
pub fn query(self: *const Tree, box: Aabb, context: anytype, comptime visit: fn (@TypeOf(context), u32) bool) void {
    if (self.root == null_node) return;
    var stack: [stack_depth]u32 = undefined;
    var top: usize = 1;
    stack[0] = self.root;
    while (top > 0) {
        top -= 1;
        const node = &self.nodes.items[stack[top]];
        if (!node.box.overlaps(box)) continue;
        if (node.isLeaf()) {
            if (!visit(context, node.value)) return;
        } else if (top + 2 <= stack_depth) {
            stack[top] = node.child1;
            stack[top + 1] = node.child2;
            top += 2;
        }
    }
}

/// Walk the leaves a ray from `origin` along `translation` may hit, calling
/// `visit(context, value, max_fraction)` for each. The visitor answers with
/// the fraction it hit that leaf at, which shortens the ray for the rest of
/// the walk; a negative number for a miss; or nought to stop.
pub fn rayCast(
    self: *const Tree,
    origin: Vec3,
    translation: Vec3,
    max_fraction: f32,
    context: anytype,
    comptime visit: fn (@TypeOf(context), u32, f32) f32,
) void {
    self.sweepCast(origin, translation, .zero, max_fraction, context, visit);
}

/// As `rayCast`, for a box of half size `half` swept along the ray: every
/// leaf box is grown by it first. What a shape cast walks.
pub fn sweepCast(
    self: *const Tree,
    origin: Vec3,
    translation: Vec3,
    half: Vec3,
    max_fraction: f32,
    context: anytype,
    comptime visit: fn (@TypeOf(context), u32, f32) f32,
) void {
    if (self.root == null_node) return;
    var fraction = max_fraction;
    var stack: [stack_depth]u32 = undefined;
    var top: usize = 1;
    stack[0] = self.root;
    while (top > 0) {
        top -= 1;
        const node = &self.nodes.items[stack[top]];
        const grown: Aabb = .{ .min = node.box.min.sub(half), .max = node.box.max.add(half) };
        if (grown.rayCast(origin, translation, fraction) == null) continue;
        if (node.isLeaf()) {
            const hit = visit(context, node.value, fraction);
            if (hit == 0) return;
            if (hit > 0 and hit < fraction) fraction = hit;
        } else if (top + 2 <= stack_depth) {
            stack[top] = node.child1;
            stack[top + 1] = node.child2;
            top += 2;
        }
    }
}

// -------------------------------------------------------------------------
// Tests
// -------------------------------------------------------------------------

fn validate(self: *const Tree) !void {
    if (self.root == null_node) return;
    try testing.expectEqual(null_node, self.nodes.items[self.root].parent);
    var leaves: u32 = 0;
    try self.validateNode(self.root, &leaves);
    try testing.expectEqual(self.leaf_count, leaves);
}

fn validateNode(self: *const Tree, index: u32, leaves: *u32) !void {
    const node = &self.nodes.items[index];
    if (node.isLeaf()) {
        try testing.expectEqual(@as(i32, 0), node.height);
        leaves.* += 1;
        return;
    }
    const a = &self.nodes.items[node.child1];
    const b = &self.nodes.items[node.child2];
    try testing.expectEqual(index, a.parent);
    try testing.expectEqual(index, b.parent);
    try testing.expectEqual(1 + @max(a.height, b.height), node.height);
    try testing.expect(node.box.contains(a.box) and node.box.contains(b.box));
    try self.validateNode(node.child1, leaves);
    try self.validateNode(node.child2, leaves);
}

const Found = struct {
    values: std.ArrayList(u32) = .empty,

    fn add(self: *Found, value: u32) bool {
        self.values.append(testing.allocator, value) catch return false;
        return true;
    }
};

fn cube(x: f32, y: f32, z: f32) Aabb {
    return .{ .min = .init(x, y, z), .max = .init(x + 1, y + 1, z + 1) };
}

test "a grid of boxes put in row by row stays balanced, and a query finds just the ones it overlaps" {
    var tree: Tree = .empty;
    defer tree.deinit(testing.allocator);
    var proxies: [512]u32 = undefined;
    for (0..8) |x| for (0..8) |y| for (0..8) |z| {
        const at = x * 64 + y * 8 + z;
        proxies[at] = try tree.insert(testing.allocator, cube(@floatFromInt(x * 2), @floatFromInt(y * 2), @floatFromInt(z * 2)), @intCast(at));
    };
    try tree.validate();
    try testing.expect(tree.height() <= 14);

    var found: Found = .{};
    defer found.values.deinit(testing.allocator);
    tree.query(.{ .min = .init(1.5, 1.5, 1.5), .max = .init(2.5, 2.5, 2.5) }, &found, Found.add);
    // Only the box at (2, 2, 2).
    try testing.expectEqual(@as(usize, 1), found.values.items.len);
    try testing.expectEqual(@as(u32, 64 + 8 + 1), found.values.items[0]);

    for (proxies[0..256]) |proxy| tree.remove(proxy);
    try tree.validate();
    try testing.expectEqual(@as(u32, 256), tree.leaf_count);
    try tree.move(testing.allocator, proxies[300], cube(100, 100, 100));
    try tree.validate();
}

const Nearest = struct {
    fn visit(_: *Nearest, value: u32, max: f32) f32 {
        // Each box is hit where its x starts, along a ray from x = -10.
        const fraction = (@as(f32, @floatFromInt(value * 2)) + 10) / 100;
        return if (fraction <= max) fraction else -1;
    }
};

test "a ray walks only what it may hit, nearer hits shortening it" {
    var tree: Tree = .empty;
    defer tree.deinit(testing.allocator);
    for (0..20) |i| _ = try tree.insert(testing.allocator, cube(@floatFromInt(i * 2), 0, 0), @intCast(i));
    var nearest: Nearest = .{};
    var best: f32 = 1;
    const Recorder = struct {
        inner: *Nearest,
        best: *f32,
        fn visit(self: @This(), value: u32, max: f32) f32 {
            const hit = self.inner.visit(value, max);
            if (hit > 0) self.best.* = @min(self.best.*, hit);
            return hit;
        }
    };
    tree.rayCast(.init(-10, 0.5, 0.5), .init(100, 0, 0), 1, Recorder{ .inner = &nearest, .best = &best }, Recorder.visit);
    try testing.expectApproxEqAbs(@as(f32, 0.1), best, 1e-6);
}

// SPDX-License-Identifier: BSD-2-Clause

//! A mesh of triangles to collide with: the level, a terrain, a room - what
//! never moves. Its triangles are kept in a tree of their own, so a body
//! near one corner of a level of a hundred thousand triangles looks at the
//! few under it.
//!
//! ```zig
//! var level = try Mesh.init(gpa, vertices, indices);
//! defer level.deinit(gpa);
//! const ground = try world.createBody(.{ .type = .static });
//! _ = try world.addShape(ground, .{ .geometry = .{ .mesh = &level } });
//! ```
//!
//! **Only on a static body.** A mesh is a surface, not a solid: it has no
//! inside, so no mass, and two of them cannot tell how far one is sunk in
//! the other. What moves is a box, a ball, a capsule, a cylinder or a hull.
//!
//! **Each triangle is solid from both sides**: what comes at it from
//! behind is stopped too, so a level whose triangles face every which way
//! still holds.

const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Aabb = geometry.Aabb;
const Tree = @import("tree.zig");

pub const Mesh = @This();

vertices: []Vec3,
triangles: [][3]u32,
tree: Tree,
bounds: Aabb,

pub const Error = Allocator.Error || error{
    /// An index past the last vertex, or a count of indices that is not a
    /// whole number of triangles.
    BadMesh,
};

/// A mesh of `vertices`, three `indices` a triangle, copied.
pub fn init(gpa: Allocator, vertices: []const Vec3, indices: []const u32) Error!Mesh {
    if (indices.len % 3 != 0) return error.BadMesh;
    for (indices) |i| if (i >= vertices.len) return error.BadMesh;
    const own_vertices = try gpa.dupe(Vec3, vertices);
    errdefer gpa.free(own_vertices);
    const triangles = try gpa.alloc([3]u32, indices.len / 3);
    errdefer gpa.free(triangles);
    var kept: usize = 0;
    for (0..triangles.len) |t| {
        const tri: [3]u32 = .{ indices[t * 3], indices[t * 3 + 1], indices[t * 3 + 2] };
        // A triangle with no area is nothing to stand on, and has no normal.
        const a = vertices[tri[0]];
        if (vertices[tri[1]].sub(a).cross(vertices[tri[2]].sub(a)).lenSq() < 1e-20) continue;
        triangles[kept] = tri;
        kept += 1;
    }
    var tree: Tree = .empty;
    errdefer tree.deinit(gpa);
    var bounds = Aabb.empty;
    for (triangles[0..kept], 0..) |tri, t| {
        const box = Aabb.ofPoints(&.{ own_vertices[tri[0]], own_vertices[tri[1]], own_vertices[tri[2]] });
        bounds = bounds.join(box);
        _ = try tree.insert(gpa, box, @intCast(t));
    }
    const shrunk = if (kept == triangles.len) triangles else blk: {
        const copy = try gpa.dupe([3]u32, triangles[0..kept]);
        gpa.free(triangles);
        break :blk copy;
    };
    return .{ .vertices = own_vertices, .triangles = shrunk, .tree = tree, .bounds = bounds };
}

pub fn deinit(self: *Mesh, gpa: Allocator) void {
    gpa.free(self.vertices);
    gpa.free(self.triangles);
    self.tree.deinit(gpa);
    self.* = undefined;
}

/// The corners of triangle `t`, in the mesh's own frame.
pub fn triangle(self: *const Mesh, t: u32) [3]Vec3 {
    const tri = self.triangles[t];
    return .{ self.vertices[tri[0]], self.vertices[tri[1]], self.vertices[tri[2]] };
}

/// Call `visit(context, t)` for each triangle whose box overlaps `box`, in
/// the mesh's own frame, until it returns false.
pub fn query(self: *const Mesh, box: Aabb, context: anytype, comptime visit: fn (@TypeOf(context), u32) bool) void {
    self.tree.query(box, context, visit);
}

pub const RayHit = struct {
    fraction: f32,
    /// Out of the side the ray came from.
    normal: Vec3,
    triangle: u32,
};

/// Where a ray, in the mesh's own frame, first meets a triangle.
pub fn rayCast(self: *const Mesh, origin: Vec3, translation: Vec3, max_fraction: f32) ?RayHit {
    const Walk = struct {
        mesh: *const Mesh,
        origin: Vec3,
        translation: Vec3,
        best: ?RayHit = null,

        fn visit(walk: *@This(), t: u32, max: f32) f32 {
            const corners = walk.mesh.triangle(t);
            const hit = rayTriangle(walk.origin, walk.translation, corners[0], corners[1], corners[2]) orelse return -1;
            if (hit.fraction > max) return -1;
            walk.best = .{ .fraction = hit.fraction, .normal = hit.normal, .triangle = t };
            return hit.fraction;
        }
    };
    var walk: Walk = .{ .mesh = self, .origin = origin, .translation = translation };
    self.tree.rayCast(origin, translation, max_fraction, &walk, Walk.visit);
    return walk.best;
}

/// Where a ray meets triangle `a b c`, from either side - Möller and
/// Trumbore's test - with the normal facing back along the ray.
pub fn rayTriangle(origin: Vec3, translation: Vec3, a: Vec3, b: Vec3, c: Vec3) ?struct { fraction: f32, normal: Vec3 } {
    const e1 = b.sub(a);
    const e2 = c.sub(a);
    const p = translation.cross(e2);
    const det = e1.dot(p);
    if (@abs(det) < 1e-12) return null;
    const inv = 1 / det;
    const s = origin.sub(a);
    const u = s.dot(p) * inv;
    if (u < 0 or u > 1) return null;
    const q = s.cross(e1);
    const v = translation.dot(q) * inv;
    if (v < 0 or u + v > 1) return null;
    const t = e2.dot(q) * inv;
    if (t < 0) return null;
    var n = e1.cross(e2).norm();
    if (n.dot(translation) > 0) n = n.neg();
    return .{ .fraction = t, .normal = n };
}

test "a mesh finds the triangles near a box and the first a ray meets, from either side" {
    // Two floors, a unit apart, each a square of two triangles.
    const vertices = [_]Vec3{
        .init(-1, 0, -1), .init(1, 0, -1), .init(1, 0, 1), .init(-1, 0, 1),
        .init(-1, 1, -1), .init(1, 1, -1), .init(1, 1, 1), .init(-1, 1, 1),
    };
    const indices = [_]u32{ 0, 2, 1, 0, 3, 2, 4, 6, 5, 4, 7, 6, 0, 0, 1 };
    var mesh = try Mesh.init(testing.allocator, &vertices, &indices);
    defer mesh.deinit(testing.allocator);
    // The triangle with no area is left out.
    try testing.expectEqual(@as(usize, 4), mesh.triangles.len);

    const Count = struct {
        n: usize = 0,
        fn add(self: *@This(), _: u32) bool {
            self.n += 1;
            return true;
        }
    };
    var near: Count = .{};
    mesh.query(.{ .min = .init(-0.1, -0.1, -0.1), .max = .init(0.1, 0.1, 0.1) }, &near, Count.add);
    try testing.expectEqual(@as(usize, 2), near.n);

    const down = mesh.rayCast(.init(0.2, 5, 0.3), .init(0, -10, 0), 1).?;
    try testing.expectApproxEqAbs(@as(f32, 0.4), down.fraction, 1e-6);
    try testing.expect(down.normal.approxEql(.init(0, 1, 0)));
    const up = mesh.rayCast(.init(0.2, -5, 0.3), .init(0, 10, 0), 1).?;
    try testing.expectApproxEqAbs(@as(f32, 0.5), up.fraction, 1e-6);
    try testing.expect(up.normal.approxEql(.init(0, -1, 0)));
    try testing.expect(mesh.rayCast(.init(3, 5, 0), .init(0, -10, 0), 1) == null);
    try testing.expectError(error.BadMesh, Mesh.init(testing.allocator, &vertices, &.{ 0, 1, 9 }));
}

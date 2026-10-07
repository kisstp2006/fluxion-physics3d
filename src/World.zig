// SPDX-License-Identifier: BSD-2-Clause

//! The bodies, the shapes on them, and one step of time.
//!
//! ```zig
//! var world: World = .init(gpa, .{});
//! defer world.deinit();
//! const ground = try world.createBody(.{ .type = .static });
//! _ = try world.addShape(ground, .box(.init(10, 0.5, 10)));
//! const crate = try world.createBody(.{ .position = .init(0, 4, 0) });
//! _ = try world.addShape(crate, .box(.init(0.5, 0.5, 0.5)));
//! try world.step(1.0 / 60.0);
//! ```
//!
//! **A step**, in order: what moves is noted and its boxes moved in the
//! tree of moving shapes; each moving shape asks both trees - the moving
//! one and the one of what never moves - what is near it; each pair near
//! enough is asked where it touches (`collide`), a mesh triangle by
//! triangle; the touching pairs become constraints (`contact`), warm from
//! what the same points pushed last step; the step is cut into substeps,
//! each integrating the velocities, solving, moving and relaxing; the
//! bounces are given; what began and stopped touching is said; and islands
//! of bodies that have been still long enough fall asleep.
//!
//! **Units are metres, kilograms and seconds**, and gravity pulls towards
//! `-y`, as the engine's 3D world is laid out.

const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const id = @import("fluxion_id");

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Quat = geometry.Quat;
const Transform = geometry.Transform;
const Aabb = geometry.Aabb;
const shape_mod = @import("shape.zig");
const Shape = shape_mod.Shape;
const Geometry = shape_mod.Geometry;
const Body = @import("body.zig");
const Tree = @import("tree.zig");
const collide = @import("collide.zig");
const Manifold = collide.Manifold;
const contact = @import("contact.zig");
const gjk = @import("gjk.zig");
const Softness = @import("Softness.zig");
const Mesh = @import("mesh.zig");

const World = @This();

pub const BodyId = Body.Id;
pub const ShapeId = Body.ShapeId;

pub const Settings = struct {
    gravity: Vec3 = .init(0, -9.81, 0),
    /// How many times a step is cut: more is stiffer and dearer.
    substeps: u32 = 4,
    /// How deep a resting contact sits: some, or it flickers on and off.
    linear_slop: f32 = 0.005,
    /// How far apart two shapes may be and still have points: what stops a
    /// fast body at a surface before it goes through.
    speculative_distance: f32 = 0.04,
    /// How far past a moving shape its box in the tree reaches, so the box
    /// need not move every step.
    aabb_margin: f32 = 0.1,
    /// The spring sinking in is taken back with. See `Softness`.
    contact_hertz: f32 = 30,
    contact_damping_ratio: f32 = 10,
    /// The fastest sinking in is taken back, in metres a second.
    max_push_speed: f32 = 3,
    /// Slower than this, a hit does not bounce.
    restitution_threshold: f32 = 1,
    /// The fastest anything goes, in metres a second.
    max_speed: f32 = 100,
    enable_sleep: bool = true,
    /// Slower than these, in metres and radians a second, a body is still.
    sleep_threshold: f32 = 0.05,
    angular_sleep_threshold: f32 = 0.05,
    /// How long an island must be still to fall asleep.
    time_to_sleep: f32 = 0.5,
    filter_rule: shape_mod.FilterRule = .both,
    friction_mix: shape_mod.Mix = .geometric_mean,
    restitution_mix: shape_mod.Mix = .maximum,
};

/// A shape on a body, as the world keeps it.
pub const ShapeEntry = struct {
    def: Shape,
    body: BodyId,
    body_index: u32,
    /// The next shape on the same body.
    next: ShapeId,
    proxy: u32,
    /// In the tree of what never moves, or of what does.
    in_static: bool,
    /// Its box in its tree: grown, for a moving shape.
    box: Aabb,
};

pub const ContactEvent = struct {
    shape_a: ShapeId,
    shape_b: ShapeId,
    body_a: BodyId,
    body_b: BodyId,
    /// One of them is a sensor: seen, not pushed.
    sensor: bool,
};

pub const RayHit = struct {
    shape: ShapeId,
    body: BodyId,
    point: Vec3,
    /// Out of the surface, back along the ray.
    normal: Vec3,
    fraction: f32,
};

/// What a ray, a cast or an overlap looks at.
pub const QueryFilter = struct {
    /// The layers it sees.
    mask: u32 = 0xFFFF_FFFF,
    /// A body passed over: the one it is cast from.
    ignore: BodyId = .none,
    /// Whether sensors are seen.
    sensors: bool = false,
    /// How far short a cast stops, and how near counts as sunk in to
    /// `penetrations`: what a character keeps clear of walls by. Nought is
    /// the slop.
    margin: f32 = 0,
};

pub const ShapeHit = struct {
    shape: ShapeId,
    body: BodyId,
    point: Vec3,
    /// Out of what it hit, towards what was cast.
    normal: Vec3,
    fraction: f32,
    /// It overlapped at the start: no normal to tell.
    initially_overlapping: bool,
};

/// How far a shape has sunk into another, and which way out.
pub const Penetration = struct {
    shape: ShapeId,
    body: BodyId,
    /// Out of the other shape, towards the one asked about.
    normal: Vec3,
    depth: f32,
    point: Vec3,
};

pub const Error = Allocator.Error;

pub const ShapeError = Allocator.Error || error{
    NoSuchBody,
    /// A mesh has no inside: it can be only on a static body.
    MeshOnMovingBody,
};

gpa: Allocator,
settings: Settings,

bodies: id.Table(Body) = .empty,
shapes: id.Table(ShapeEntry) = .empty,
static_tree: Tree = .empty,
moving_tree: Tree = .empty,
/// Pairs of bodies that never collide, by their handles.
exceptions: std.AutoHashMapUnmanaged([2]u64, void) = .empty,

/// What touches now, and what touched at the end of the last step.
contacts: ContactMap = .empty,
previous: ContactMap = .empty,
begin_events: std.ArrayList(ContactEvent) = .empty,
end_events: std.ArrayList(ContactEvent) = .empty,
/// The touches of shapes taken out since the last step: they end in the
/// next one.
dropped_ends: std.ArrayList(ContactEvent) = .empty,

// Worked out each step, kept for their memory.
movers: std.ArrayList(u32) = .empty,
mover_mark: std.ArrayList(u64) = .empty,
tight: std.ArrayList(Aabb) = .empty,
pairs: std.ArrayList([2]u32) = .empty,
found: std.ArrayList(Found) = .empty,
constraints: std.ArrayList(contact.Constraint) = .empty,
island_parent: std.ArrayList(u32) = .empty,
island_rest: std.ArrayList(f32) = .empty,
step_count: u64 = 0,

const ContactMap = std.AutoHashMapUnmanaged(contact.PairKey, Stored);

/// What is kept of a pair that touches.
const Stored = struct {
    shape_a: ShapeId,
    shape_b: ShapeId,
    body_a: BodyId,
    body_b: BodyId,
    sensor: bool,
    impulses: [collide.max_points]contact.Impulse = @splat(.{}),
};

/// A pair's manifold, found this step.
const Found = struct {
    key: contact.PairKey,
    a: u32,
    b: u32,
    manifold: Manifold,
    sensor: bool,
};

pub fn init(gpa: Allocator, settings: Settings) World {
    return .{ .gpa = gpa, .settings = settings };
}

pub fn deinit(self: *World) void {
    const gpa = self.gpa;
    self.bodies.deinit(gpa);
    self.shapes.deinit(gpa);
    self.static_tree.deinit(gpa);
    self.moving_tree.deinit(gpa);
    self.exceptions.deinit(gpa);
    self.contacts.deinit(gpa);
    self.previous.deinit(gpa);
    self.begin_events.deinit(gpa);
    self.end_events.deinit(gpa);
    self.dropped_ends.deinit(gpa);
    self.movers.deinit(gpa);
    self.mover_mark.deinit(gpa);
    self.tight.deinit(gpa);
    self.pairs.deinit(gpa);
    self.found.deinit(gpa);
    self.constraints.deinit(gpa);
    self.island_parent.deinit(gpa);
    self.island_rest.deinit(gpa);
    self.* = undefined;
}

// -------------------------------------------------------------------------
// Bodies and shapes
// -------------------------------------------------------------------------

pub fn createBody(self: *World, def: Body.Def) Error!BodyId {
    var made = Body.init(def);
    if (made.type == .dynamic) {
        // Until its shapes say otherwise: a unit ball's mass.
        made.mass = 1;
        made.inv_mass = 1;
        made.inertia = geometry.diagonal(.splat(0.4));
        made.inv_inertia_local = geometry.diagonal(.splat(2.5));
    }
    made.updateInertia();
    return self.bodies.add(self.gpa, made);
}

pub fn body(self: *World, handle: BodyId) ?*Body {
    return self.bodies.get(handle);
}

pub fn bodyConst(self: *const World, handle: BodyId) ?*const Body {
    return self.bodies.getConst(handle);
}

pub fn bodyCount(self: *const World) usize {
    return self.bodies.count();
}

fn bodyAt(self: *World, index: u32) *Body {
    return &self.bodies.slots.items[index].value.?;
}

fn handleOfBody(self: *const World, index: u32) BodyId {
    return .{ .index = index, .generation = self.bodies.slots.items[index].generation };
}

fn handleOfShape(self: *const World, index: u32) ShapeId {
    return .{ .index = index, .generation = self.shapes.slots.items[index].generation };
}

fn shapeAt(self: *World, index: u32) *ShapeEntry {
    return &self.shapes.slots.items[index].value.?;
}

/// Take a body away, with its shapes. What it touched stops touching it.
pub fn destroyBody(self: *World, handle: BodyId) void {
    const b = self.bodies.get(handle) orelse return;
    var cursor = b.first_shape;
    while (self.shapes.get(cursor)) |entry| {
        const next = entry.next;
        self.dropShape(cursor, entry);
        cursor = next;
    }
    _ = self.bodies.remove(handle);
}

/// Put a shape on a body. A body's mass is worked out again from its
/// shapes, and it wakes.
pub fn addShape(self: *World, body_handle: BodyId, def: Shape) ShapeError!ShapeId {
    const b = self.bodies.get(body_handle) orelse return error.NoSuchBody;
    if (def.geometry == .mesh and b.type != .static) return error.MeshOnMovingBody;
    const xf = b.transform.mul(def.offset);
    const in_static = b.type == .static;
    var box = def.geometry.aabb(xf);
    if (!in_static) box = box.grow(self.settings.aabb_margin);
    const handle = try self.shapes.add(self.gpa, .{
        .def = def,
        .body = body_handle,
        .body_index = body_handle.index,
        .next = b.first_shape,
        .proxy = 0,
        .in_static = in_static,
        .box = box,
    });
    const tree = if (in_static) &self.static_tree else &self.moving_tree;
    const proxy = tree.insert(self.gpa, box, handle.index) catch |err| {
        _ = self.shapes.remove(handle);
        return err;
    };
    self.shapes.get(handle).?.proxy = proxy;
    const owner = self.bodies.get(body_handle).?;
    owner.first_shape = handle;
    owner.shape_count += 1;
    self.updateMass(body_handle);
    owner.wake();
    return handle;
}

pub fn removeShape(self: *World, handle: ShapeId) void {
    const entry = self.shapes.get(handle) orelse return;
    const owner = entry.body;
    // Out of its body's list.
    if (self.bodies.get(owner)) |b| {
        if (b.first_shape.eql(handle)) {
            b.first_shape = entry.next;
        } else {
            var cursor = b.first_shape;
            while (self.shapes.get(cursor)) |e| : (cursor = e.next) {
                if (e.next.eql(handle)) {
                    e.next = entry.next;
                    break;
                }
            }
        }
        b.shape_count -= 1;
    }
    self.dropShape(handle, entry);
    self.updateMass(owner);
    if (self.bodies.get(owner)) |b| b.wake();
}

/// Out of its tree and its contacts, and gone. What it touched ends in the
/// next step, and what it held up or leaned on wakes.
fn dropShape(self: *World, handle: ShapeId, entry: *ShapeEntry) void {
    const tree = if (entry.in_static) &self.static_tree else &self.moving_tree;
    tree.remove(entry.proxy);
    for ([_]*ContactMap{ &self.contacts, &self.previous }) |map| {
        while (true) {
            var it = map.iterator();
            const gone = while (it.next()) |kv| {
                if (kv.value_ptr.shape_a.eql(handle) or kv.value_ptr.shape_b.eql(handle)) break kv;
            } else null;
            const kv = gone orelse break;
            if (map == &self.contacts) self.dropTouch(kv.value_ptr.*);
            _ = map.remove(kv.key_ptr.*);
        }
    }
    _ = self.shapes.remove(handle);
}

fn dropTouch(self: *World, stored: Stored) void {
    if (!stored.sensor) {
        if (self.bodies.get(stored.body_a)) |b| b.wake();
        if (self.bodies.get(stored.body_b)) |b| b.wake();
    }
    // A mesh's triangles are one touch.
    for (self.dropped_ends.items) |e| if (e.shape_a.eql(stored.shape_a) and e.shape_b.eql(stored.shape_b)) return;
    // Out of memory, the end goes unsaid: the shape is gone either way.
    self.dropped_ends.append(self.gpa, eventOf(stored)) catch {};
}

pub fn shape(self: *World, handle: ShapeId) ?*ShapeEntry {
    return self.shapes.get(handle);
}

pub fn shapeCount(self: *const World) usize {
    return self.shapes.count();
}

/// Where a shape is in the world.
pub fn shapeTransform(self: *World, entry: *const ShapeEntry) Transform {
    return self.bodyAt(entry.body_index).transform.mul(entry.def.offset);
}

/// Move a body there, at once: it wakes, and what it now overlaps is
/// pushed out of it by the next step.
pub fn setTransform(self: *World, handle: BodyId, position: Vec3, rotation: Quat) void {
    const b = self.bodies.get(handle) orelse return;
    b.transform = .{ .position = position, .rotation = rotation.norm() };
    b.center = b.transform.apply(b.local_center);
    b.updateInertia();
    b.wake();
    var cursor = b.first_shape;
    while (self.shapes.get(cursor)) |entry| : (cursor = entry.next) {
        var box = entry.def.geometry.aabb(b.transform.mul(entry.def.offset));
        if (!entry.in_static) box = box.grow(self.settings.aabb_margin);
        entry.box = box;
        const tree = if (entry.in_static) &self.static_tree else &self.moving_tree;
        tree.move(self.gpa, entry.proxy, box) catch {};
    }
}

/// Work a body's mass, centre and inertia out again from its shapes'
/// densities.
pub fn updateMass(self: *World, handle: BodyId) void {
    const b = self.bodies.get(handle) orelse return;
    b.mass = 0;
    b.inv_mass = 0;
    b.inertia = .zero;
    b.inv_inertia_local = .zero;
    b.local_center = .zero;
    if (b.type != .dynamic) {
        b.center = b.transform.position;
        b.updateInertia();
        return;
    }
    // The centre first, then each shape's inertia moved to it.
    var total: f32 = 0;
    var moment: Vec3 = .zero;
    var cursor = b.first_shape;
    while (self.shapes.get(cursor)) |entry| : (cursor = entry.next) {
        if (entry.def.sensor) continue;
        const m = entry.def.geometry.massOf(entry.def.material.density);
        total += m.mass;
        moment = moment.mulAdd(entry.def.offset.apply(m.center), m.mass);
    }
    if (total <= 0) {
        b.mass = 1;
        b.inv_mass = 1;
        b.inertia = geometry.diagonal(.splat(0.4));
    } else {
        const center = moment.scale(1 / total);
        var inertia: geometry.Mat3 = .zero;
        cursor = b.first_shape;
        while (self.shapes.get(cursor)) |entry| : (cursor = entry.next) {
            if (entry.def.sensor) continue;
            const m = entry.def.geometry.massOf(entry.def.material.density);
            const turned = geometry.rotateTensor(geometry.rotationOf(entry.def.offset.rotation), m.inertia);
            const d = entry.def.offset.apply(m.center).sub(center);
            const shift = geometry.addMat(geometry.diagonal(.splat(d.dot(d))), geometry.scaleMat(geometry.outer(d, d), -1));
            inertia = geometry.addMat(inertia, geometry.addMat(turned, geometry.scaleMat(shift, m.mass)));
        }
        b.mass = total;
        b.inv_mass = 1 / total;
        b.inertia = inertia;
        b.local_center = center;
    }
    b.inv_inertia_local = if (b.lock_rotation) .zero else geometry.inverseOrZero(b.inertia);
    b.center = b.transform.apply(b.local_center);
    b.updateInertia();
}

/// Two bodies that never collide, whatever their filters say.
pub fn addCollisionException(self: *World, a: BodyId, b: BodyId) Error!void {
    try self.exceptions.put(self.gpa, exceptionKey(a, b), {});
}

pub fn removeCollisionException(self: *World, a: BodyId, b: BodyId) void {
    _ = self.exceptions.remove(exceptionKey(a, b));
}

fn exceptionKey(a: BodyId, b: BodyId) [2]u64 {
    const x = a.toInt();
    const y = b.toInt();
    return if (x < y) .{ x, y } else .{ y, x };
}

// -------------------------------------------------------------------------
// The step
// -------------------------------------------------------------------------

pub fn step(self: *World, dt: f32) Error!void {
    if (dt <= 0) return;
    const gpa = self.gpa;
    const s = self.settings;
    self.step_count += 1;
    self.begin_events.clearRetainingCapacity();
    self.end_events.clearRetainingCapacity();
    try self.end_events.appendSlice(gpa, self.dropped_ends.items);
    self.dropped_ends.clearRetainingCapacity();

    const substeps = @max(s.substeps, 1);
    const h = dt / @as(f32, @floatFromInt(substeps));
    const inv_h = 1 / h;
    const hertz = @min(s.contact_hertz, 0.25 * inv_h);
    const softness: Softness = .of(hertz, s.contact_damping_ratio, h);
    const static_softness: Softness = .of(2 * hertz, s.contact_damping_ratio, h);
    const solver_step: contact.Step = .{ .inv_h = inv_h, .slop = s.linear_slop, .max_push = s.max_push_speed };

    // 1. What moves this step.
    const body_slots = self.bodies.slotCount();
    try self.mover_mark.resize(gpa, body_slots);
    self.movers.clearRetainingCapacity();
    for (self.bodies.slots.items, 0..) |*slot, i| {
        const b = &(slot.value orelse continue);
        if (b.type == .dynamic and !b.awake and (!b.allow_sleep or !s.enable_sleep)) b.wake();
        if (!b.isActive()) continue;
        b.updateInertia();
        try self.markMover(@intCast(i));
    }

    // 2. Their shapes' boxes, swept along where they are going, and moved in
    //    the tree when they have left the grown box it has.
    try self.tight.resize(gpa, self.shapes.slotCount());
    for (self.movers.items) |index| {
        const b = self.bodyAt(index);
        var cursor = b.first_shape;
        while (self.shapes.get(cursor)) |entry| : (cursor = entry.next) {
            const xf = b.transform.mul(entry.def.offset);
            const box = entry.def.geometry.aabb(xf).grow(s.speculative_distance).sweep(b.linear_velocity.scale(dt));
            self.tight.items[cursor.index] = box;
            if (!entry.box.contains(box)) {
                entry.box = box.grow(s.aabb_margin);
                try self.moving_tree.move(gpa, entry.proxy, entry.box);
            }
        }
    }

    // 3. Which shapes are near which.
    self.pairs.clearRetainingCapacity();
    for (self.movers.items) |index| {
        const b = self.bodyAt(index);
        var cursor = b.first_shape;
        while (self.shapes.get(cursor)) |entry| : (cursor = entry.next) {
            var finder: PairFinder = .{ .world = self, .shape = cursor.index, .entry = entry };
            self.moving_tree.query(self.tight.items[cursor.index], &finder, PairFinder.visit);
            self.static_tree.query(self.tight.items[cursor.index], &finder, PairFinder.visit);
            if (finder.failed) return error.OutOfMemory;
        }
    }

    // 4. Where each pair touches.
    self.found.clearRetainingCapacity();
    for (self.pairs.items) |pair| try self.narrowPhase(pair[0], pair[1]);

    // 5. A sleeper something runs into wakes, and moves this step.
    for (self.found.items) |f| {
        if (f.sensor) continue;
        for ([_]u32{ f.a, f.b }) |shape_index| {
            const bi = self.shapeAt(shape_index).body_index;
            const b = self.bodyAt(bi);
            if (b.type == .dynamic and !b.awake) {
                b.wake();
                b.updateInertia();
                try self.markMover(bi);
            }
        }
    }

    // 6. The touching pairs, remembered, and their constraints.
    std.mem.swap(ContactMap, &self.contacts, &self.previous);
    self.contacts.clearRetainingCapacity();
    self.constraints.clearRetainingCapacity();
    for (self.found.items) |f| {
        const ea = self.shapeAt(f.a);
        const eb = self.shapeAt(f.b);
        var stored: Stored = .{
            .shape_a = self.handleOfShape(f.a),
            .shape_b = self.handleOfShape(f.b),
            .body_a = self.handleOfBody(ea.body_index),
            .body_b = self.handleOfBody(eb.body_index),
            .sensor = f.sensor,
        };
        if (self.previous.get(f.key)) |last| stored.impulses = last.impulses;
        if (!f.sensor) {
            const friction = s.friction_mix.of(ea.def.material.friction, eb.def.material.friction);
            const restitution = s.restitution_mix.of(ea.def.material.restitution, eb.def.material.restitution);
            const warm: ?[collide.max_points]contact.Impulse = if (self.previous.contains(f.key)) stored.impulses else null;
            try self.constraints.append(gpa, contact.prepare(f.key, &f.manifold, self.bodyAt(ea.body_index), self.bodyAt(eb.body_index), ea.body_index, eb.body_index, friction, restitution, warm, softness, static_softness));
        }
        if (self.touches(f)) try self.contacts.put(gpa, f.key, stored);
    }

    // 7. The substeps.
    for (0..substeps) |_| {
        self.integrateVelocities(h);
        for (self.constraints.items) |*c| contact.warmStart(c, self.bodyAt(c.body_a), self.bodyAt(c.body_b));
        for (self.constraints.items) |*c| contact.solve(c, self.bodyAt(c.body_a), self.bodyAt(c.body_b), solver_step, .solve);
        self.integratePositions(h);
        for (self.constraints.items) |*c| contact.solve(c, self.bodyAt(c.body_a), self.bodyAt(c.body_b), solver_step, .relax);
    }

    // 8. The bounces; and the forces, which pushed for the whole step, gone.
    for (self.constraints.items) |*c| contact.restitute(c, self.bodyAt(c.body_a), self.bodyAt(c.body_b), s.restitution_threshold);
    for (self.movers.items) |index| {
        const b = self.bodyAt(index);
        b.force = .zero;
        b.torque = .zero;
    }

    // 9. What each point pushed, for the next step.
    for (self.constraints.items) |*c| {
        if (self.contacts.getPtr(c.key)) |kept| kept.impulses = contact.impulsesOf(c);
    }

    // 10. What began and stopped touching. A pair of which nothing moved
    //     this step - both asleep, or one asleep on what never moves - was
    //     not looked at, and touches as it did.
    var it = self.previous.iterator();
    while (it.next()) |kv| {
        if (self.contacts.contains(kv.key_ptr.*)) continue;
        const was = kv.value_ptr.*;
        const still = self.shapes.contains(was.shape_a) and self.shapes.contains(was.shape_b);
        if (still and !self.isMover(was.body_a) and !self.isMover(was.body_b)) {
            try self.contacts.put(gpa, kv.key_ptr.*, was);
            continue;
        }
        if (still) try self.endOnce(was);
    }
    var now = self.contacts.iterator();
    while (now.next()) |kv| {
        if (!self.previous.contains(kv.key_ptr.*)) try self.beginOnce(kv.value_ptr.*);
    }

    // 11. Islands that have been still long enough fall asleep.
    if (s.enable_sleep) try self.updateSleep(dt);
}

fn markMover(self: *World, index: u32) Error!void {
    if (self.mover_mark.items.len <= index) try self.mover_mark.resize(self.gpa, index + 1);
    if (self.mover_mark.items[index] == self.step_count) return;
    self.mover_mark.items[index] = self.step_count;
    try self.movers.append(self.gpa, index);
}

fn isMover(self: *const World, handle: BodyId) bool {
    return handle.index < self.mover_mark.items.len and self.mover_mark.items[handle.index] == self.step_count;
}

fn eventOf(stored: Stored) ContactEvent {
    return .{ .shape_a = stored.shape_a, .shape_b = stored.shape_b, .body_a = stored.body_a, .body_b = stored.body_b, .sensor = stored.sensor };
}

/// A mesh's triangles are pairs of their own: two of them touching is one
/// touch, said once.
fn beginOnce(self: *World, stored: Stored) Error!void {
    for (self.begin_events.items) |e| if (e.shape_a.eql(stored.shape_a) and e.shape_b.eql(stored.shape_b)) return;
    if (self.wasTouching(stored, &self.previous)) return;
    try self.begin_events.append(self.gpa, eventOf(stored));
}

fn endOnce(self: *World, stored: Stored) Error!void {
    for (self.end_events.items) |e| if (e.shape_a.eql(stored.shape_a) and e.shape_b.eql(stored.shape_b)) return;
    if (self.wasTouching(stored, &self.contacts)) return;
    try self.end_events.append(self.gpa, eventOf(stored));
}

/// Whether `map` has the same two shapes touching, by any triangle.
fn wasTouching(_: *World, stored: Stored, map: *const ContactMap) bool {
    var it = map.iterator();
    while (it.next()) |kv| {
        if (kv.value_ptr.shape_a.eql(stored.shape_a) and kv.value_ptr.shape_b.eql(stored.shape_b)) return true;
    }
    return false;
}

/// Touching, not only near: sunk in, or within twice the slop of it.
fn touches(self: *const World, f: Found) bool {
    if (f.manifold.count == 0) return false;
    return f.manifold.deepest() <= if (f.sensor) 0 else 2 * self.settings.linear_slop;
}

const PairFinder = struct {
    world: *World,
    shape: u32,
    entry: *ShapeEntry,
    failed: bool = false,

    fn visit(self: *PairFinder, other: u32) bool {
        if (other == self.shape) return true;
        const world = self.world;
        const theirs = world.shapeAt(other);
        if (theirs.body_index == self.entry.body_index) return true;
        // Both moving: found from both, and kept from the lower.
        const other_moves = world.mover_mark.items.len > theirs.body_index and world.mover_mark.items[theirs.body_index] == world.step_count;
        if (other_moves and other < self.shape) return true;
        if (!world.settings.filter_rule.touches(self.entry.def.filter, theirs.def.filter)) return true;
        const mine = world.bodyAt(self.entry.body_index);
        const their_body = world.bodyAt(theirs.body_index);
        const sensor = self.entry.def.sensor or theirs.def.sensor;
        // Two that cannot both move never push each other.
        if (!sensor and mine.type != .dynamic and their_body.type != .dynamic) return true;
        if (world.exceptions.contains(exceptionKey(world.handleOfBody(self.entry.body_index), world.handleOfBody(theirs.body_index)))) return true;
        const a = @min(self.shape, other);
        const b = @max(self.shape, other);
        world.pairs.append(world.gpa, .{ a, b }) catch {
            self.failed = true;
            return false;
        };
        return true;
    }
};

fn narrowPhase(self: *World, a: u32, b: u32) Error!void {
    const ea = self.shapeAt(a);
    const eb = self.shapeAt(b);
    const xa = self.shapeTransform(ea);
    const xb = self.shapeTransform(eb);
    const sensor = ea.def.sensor or eb.def.sensor;
    const margin = if (sensor) 0 else self.settings.speculative_distance;
    const key_a = self.handleOfShape(a).toInt();
    const key_b = self.handleOfShape(b).toInt();
    if (ea.def.geometry == .mesh or eb.def.geometry == .mesh) {
        if (ea.def.geometry == .mesh and eb.def.geometry == .mesh) return;
        const mesh_first = ea.def.geometry == .mesh;
        const mesh_entry = if (mesh_first) ea else eb;
        const other_entry = if (mesh_first) eb else ea;
        const mesh_xf = if (mesh_first) xa else xb;
        const other_xf = if (mesh_first) xb else xa;
        const mesh = mesh_entry.def.geometry.mesh;
        var store_other: collide.Storage = .{};
        const other_solid = collide.solidOf(other_entry.def.geometry, other_xf, &store_other).?;
        const local_box = boxInFrame(other_entry.def.geometry.aabb(other_xf).grow(margin), mesh_xf);
        var triangles: TriangleList = .{ .gpa = self.gpa };
        defer triangles.list.deinit(self.gpa);
        mesh.query(local_box, &triangles, TriangleList.add);
        if (triangles.failed) return error.OutOfMemory;
        for (triangles.list.items) |t| {
            const corners = mesh.triangle(t);
            var store_tri: collide.Storage = .{};
            const tri = collide.triangleSolid(&store_tri, mesh_xf.apply(corners[0]), mesh_xf.apply(corners[1]), mesh_xf.apply(corners[2]));
            const m = if (mesh_first) collide.collide(tri, other_solid, margin) else collide.collide(other_solid, tri, margin);
            if (m.count == 0) continue;
            try self.found.append(self.gpa, .{ .key = .{ .a = key_a, .b = key_b, .part = t + 1 }, .a = a, .b = b, .manifold = m, .sensor = sensor });
        }
        return;
    }
    var store_a: collide.Storage = .{};
    var store_b: collide.Storage = .{};
    const sa = collide.solidOf(ea.def.geometry, xa, &store_a).?;
    const sb = collide.solidOf(eb.def.geometry, xb, &store_b).?;
    const m = collide.collide(sa, sb, margin);
    if (m.count == 0) return;
    try self.found.append(self.gpa, .{ .key = .{ .a = key_a, .b = key_b }, .a = a, .b = b, .manifold = m, .sensor = sensor });
}

const TriangleList = struct {
    gpa: Allocator,
    list: std.ArrayList(u32) = .empty,
    failed: bool = false,

    fn add(self: *TriangleList, t: u32) bool {
        self.list.append(self.gpa, t) catch {
            self.failed = true;
            return false;
        };
        return true;
    }
};

/// A world box seen from the frame `xf`: the box round its corners there.
fn boxInFrame(box: Aabb, xf: Transform) Aabb {
    var out = Aabb.empty;
    for (0..8) |i| {
        const corner = Vec3.init(
            if (i & 1 != 0) box.max.x else box.min.x,
            if (i & 2 != 0) box.max.y else box.min.y,
            if (i & 4 != 0) box.max.z else box.min.z,
        );
        out = out.add(xf.unapply(corner));
    }
    return out;
}

/// One substep's worth of gravity, forces and damping.
fn integrateVelocities(self: *World, h: f32) void {
    const s = self.settings;
    for (self.movers.items) |index| {
        const b = self.bodyAt(index);
        if (b.type != .dynamic) continue;
        const acceleration = s.gravity.scale(b.gravity_scale).mulAdd(b.force, b.inv_mass);
        b.linear_velocity = b.linear_velocity.mulAdd(acceleration, h);
        b.angular_velocity = b.angular_velocity.add(geometry.mulVec(b.inv_inertia, b.torque).scale(h));
        // Damping as `v /= 1 + c h`: stable for any `c` and any `h`.
        b.linear_velocity = b.linear_velocity.scale(1 / (1 + h * b.linear_damping));
        b.angular_velocity = b.angular_velocity.scale(1 / (1 + h * b.angular_damping));
        // Past a quarter turn a substep the solver's straight-line picture
        // of a lever arm is no use; and nothing goes faster than the most.
        const turn = h * b.angular_velocity.len();
        if (turn > 0.25 * std.math.pi) b.angular_velocity = b.angular_velocity.scale(0.25 * std.math.pi / turn);
        const speed_sq = b.linear_velocity.lenSq();
        if (speed_sq > s.max_speed * s.max_speed) b.linear_velocity = b.linear_velocity.scale(s.max_speed / @sqrt(speed_sq));
        if (b.lock_rotation) b.angular_velocity = .zero;
    }
}

fn integratePositions(self: *World, h: f32) void {
    for (self.movers.items) |index| {
        const b = self.bodyAt(index);
        if (b.type == .static) continue;
        b.center = b.center.mulAdd(b.linear_velocity, h);
        const w = b.angular_velocity;
        if (!w.eql(.zero)) {
            const q = b.transform.rotation;
            const spin = (Quat{ .x = w.x, .y = w.y, .z = w.z, .w = 0 }).mul(q);
            const half = 0.5 * h;
            b.transform.rotation = (Quat{ .x = q.x + spin.x * half, .y = q.y + spin.y * half, .z = q.z + spin.z * half, .w = q.w + spin.w * half }).norm();
            b.updateInertia();
        }
        b.placeFromCenter();
    }
}

fn findRoot(self: *World, index: u32) u32 {
    var at = index;
    while (self.island_parent.items[at] != at) {
        self.island_parent.items[at] = self.island_parent.items[self.island_parent.items[at]];
        at = self.island_parent.items[at];
    }
    return at;
}

fn updateSleep(self: *World, dt: f32) Error!void {
    const s = self.settings;
    const slots = self.bodies.slotCount();
    try self.island_parent.resize(self.gpa, slots);
    try self.island_rest.resize(self.gpa, slots);
    for (self.movers.items) |index| {
        self.island_parent.items[index] = index;
        const b = self.bodyAt(index);
        if (b.type != .dynamic) continue;
        const still = b.linear_velocity.lenSq() < s.sleep_threshold * s.sleep_threshold and
            b.angular_velocity.lenSq() < s.angular_sleep_threshold * s.angular_sleep_threshold;
        b.sleep_time = if (still and b.allow_sleep) b.sleep_time + dt else 0;
    }
    // Islands: bodies joined by what they touch. Riding something that
    // moves, a body is never still.
    for (self.constraints.items) |c| {
        const a = self.bodyAt(c.body_a);
        const b = self.bodyAt(c.body_b);
        if (a.type == .kinematic and a.isActive()) b.sleep_time = 0;
        if (b.type == .kinematic and b.isActive()) a.sleep_time = 0;
        if (a.type != .dynamic or b.type != .dynamic) continue;
        if (!self.isMover(self.handleOfBody(c.body_a)) or !self.isMover(self.handleOfBody(c.body_b))) continue;
        const ra = self.findRoot(c.body_a);
        const rb = self.findRoot(c.body_b);
        if (ra != rb) self.island_parent.items[ra] = rb;
    }
    for (self.movers.items) |index| self.island_rest.items[index] = std.math.inf(f32);
    for (self.movers.items) |index| {
        const b = self.bodyAt(index);
        if (b.type != .dynamic) continue;
        const root = self.findRoot(index);
        self.island_rest.items[root] = @min(self.island_rest.items[root], b.sleep_time);
    }
    for (self.movers.items) |index| {
        const b = self.bodyAt(index);
        if (b.type != .dynamic) continue;
        if (self.island_rest.items[self.findRoot(index)] < s.time_to_sleep) continue;
        b.awake = false;
        b.linear_velocity = .zero;
        b.angular_velocity = .zero;
    }
}

pub fn beginEvents(self: *const World) []const ContactEvent {
    return self.begin_events.items;
}

pub fn endEvents(self: *const World) []const ContactEvent {
    return self.end_events.items;
}

/// How many pairs touch now, sensors included.
pub fn touchingCount(self: *const World) usize {
    return self.contacts.count();
}

pub fn awakeCount(self: *World) usize {
    var n: usize = 0;
    var it = self.bodies.iterator();
    while (it.next()) |e| {
        if (e.value.type == .dynamic and e.value.awake) n += 1;
    }
    return n;
}

// -------------------------------------------------------------------------
// Asking
// -------------------------------------------------------------------------

fn sees(self: *World, entry: *const ShapeEntry, filter: QueryFilter) bool {
    if (entry.def.filter.layer & filter.mask == 0) return false;
    if (entry.def.sensor and !filter.sensors) return false;
    if (!filter.ignore.isNone() and filter.ignore.index == entry.body_index and self.bodies.slots.items[entry.body_index].generation == filter.ignore.generation) return false;
    return true;
}

/// The first shape a ray from `origin` along `translation` meets, from
/// outside it.
pub fn castRay(self: *World, origin: Vec3, translation: Vec3, filter: QueryFilter) ?RayHit {
    const Walk = struct {
        world: *World,
        origin: Vec3,
        translation: Vec3,
        filter: QueryFilter,
        best: ?RayHit = null,

        fn visit(walk: *@This(), index: u32, max: f32) f32 {
            const entry = walk.world.shapeAt(index);
            if (!walk.world.sees(entry, walk.filter)) return -1;
            const hit = entry.def.geometry.rayCast(walk.world.shapeTransform(entry), walk.origin, walk.translation, max) orelse return -1;
            walk.best = .{
                .shape = walk.world.handleOfShape(index),
                .body = entry.body,
                .point = walk.origin.mulAdd(walk.translation, hit.fraction),
                .normal = hit.normal,
                .fraction = hit.fraction,
            };
            return hit.fraction;
        }
    };
    var walk: Walk = .{ .world = self, .origin = origin, .translation = translation, .filter = filter };
    self.static_tree.rayCast(origin, translation, 1, &walk, Walk.visit);
    const limit = if (walk.best) |b| b.fraction else 1;
    self.moving_tree.rayCast(origin, translation, limit, &walk, Walk.visit);
    return walk.best;
}

/// Every shape a ray from `origin` along `translation` meets, from outside
/// each, nearest first - the nearest as many as `found` holds: what a
/// pointer that may pass through some things and stop at others walks.
pub fn castRayAll(self: *World, origin: Vec3, translation: Vec3, filter: QueryFilter, found: []RayHit) []RayHit {
    const Walk = struct {
        world: *World,
        origin: Vec3,
        translation: Vec3,
        filter: QueryFilter,
        found: []RayHit,
        count: usize = 0,

        fn visit(walk: *@This(), index: u32, _: f32) f32 {
            const entry = walk.world.shapeAt(index);
            if (!walk.world.sees(entry, walk.filter)) return -1;
            const hit = entry.def.geometry.rayCast(walk.world.shapeTransform(entry), walk.origin, walk.translation, 1) orelse return -1;
            const met: RayHit = .{
                .shape = walk.world.handleOfShape(index),
                .body = entry.body,
                .point = walk.origin.mulAdd(walk.translation, hit.fraction),
                .normal = hit.normal,
                .fraction = hit.fraction,
            };
            if (walk.count < walk.found.len) {
                walk.found[walk.count] = met;
                walk.count += 1;
            } else if (walk.found.len > 0) {
                // Full: the nearest are kept, in place of the furthest.
                var furthest: usize = 0;
                for (walk.found, 0..) |kept, i| {
                    if (kept.fraction > walk.found[furthest].fraction) furthest = i;
                }
                if (met.fraction < walk.found[furthest].fraction) walk.found[furthest] = met;
            }
            // Not shortened: every one along the whole ray.
            return -1;
        }
    };
    var walk: Walk = .{ .world = self, .origin = origin, .translation = translation, .filter = filter, .found = found };
    self.static_tree.rayCast(origin, translation, 1, &walk, Walk.visit);
    self.moving_tree.rayCast(origin, translation, 1, &walk, Walk.visit);
    const hits = found[0..walk.count];
    std.mem.sort(RayHit, hits, {}, nearer);
    return hits;
}

fn nearer(_: void, a: RayHit, b: RayHit) bool {
    return a.fraction < b.fraction;
}

/// Where `g`, placed by `xf` and moved along `translation`, first comes to
/// touch a shape - stopping a little short, the slop away.
pub fn castShape(self: *World, g: Geometry, xf: Transform, translation: Vec3, filter: QueryFilter) ?ShapeHit {
    var store: collide.Storage = .{};
    const moving = collide.proxyOf(g, xf, &store) orelse return null;
    const box = g.aabb(xf);
    const Walk = struct {
        world: *World,
        moving: gjk.Proxy,
        translation: Vec3,
        filter: QueryFilter,
        swept: Aabb,
        best: ?ShapeHit = null,

        fn visit(walk: *@This(), index: u32, max: f32) f32 {
            const world = walk.world;
            const entry = world.shapeAt(index);
            if (!world.sees(entry, walk.filter)) return -1;
            const shape_xf = world.shapeTransform(entry);
            const target = @max(walk.filter.margin, world.settings.linear_slop);
            var best: ?gjk.CastResult = null;
            if (entry.def.geometry == .mesh) {
                const mesh = entry.def.geometry.mesh;
                var triangles: TriangleList = .{ .gpa = world.gpa };
                defer triangles.list.deinit(world.gpa);
                mesh.query(boxInFrame(walk.swept, shape_xf), &triangles, TriangleList.add);
                for (triangles.list.items) |t| {
                    const corners = mesh.triangle(t);
                    var store_tri: collide.Storage = .{};
                    const tri = collide.triangleProxy(&store_tri, shape_xf.apply(corners[0]), shape_xf.apply(corners[1]), shape_xf.apply(corners[2]));
                    const hit = gjk.cast(tri, walk.moving, walk.translation, target, max) orelse continue;
                    if (best == null or hit.fraction < best.?.fraction) best = hit;
                }
            } else {
                var other_store: collide.Storage = .{};
                const other = collide.proxyOf(entry.def.geometry, shape_xf, &other_store).?;
                best = gjk.cast(other, walk.moving, walk.translation, target, max);
            }
            const hit = best orelse return -1;
            if (walk.best != null and hit.fraction >= walk.best.?.fraction) return -1;
            walk.best = .{
                .shape = world.handleOfShape(index),
                .body = entry.body,
                .point = hit.point,
                .normal = hit.normal,
                .fraction = hit.fraction,
                .initially_overlapping = hit.initially_overlapping,
            };
            return @max(hit.fraction, 1e-6);
        }
    };
    var walk: Walk = .{ .world = self, .moving = moving, .translation = translation, .filter = filter, .swept = box.sweep(translation) };
    const half = box.extent();
    self.static_tree.sweepCast(box.center(), translation, half, 1, &walk, Walk.visit);
    const limit = if (walk.best) |b| @max(b.fraction, 1e-6) else 1;
    self.moving_tree.sweepCast(box.center(), translation, half, limit, &walk, Walk.visit);
    return walk.best;
}

/// The shapes `g`, placed by `xf`, overlaps, into `found`.
pub fn overlapShape(self: *World, g: Geometry, xf: Transform, filter: QueryFilter, found: []ShapeId) []ShapeId {
    var deep: [32]Penetration = undefined;
    var count: usize = 0;
    for (self.penetrations(g, xf, filter, &deep)) |p| {
        const seen = for (found[0..count]) |f| {
            if (f.eql(p.shape)) break true;
        } else false;
        if (seen or count == found.len) continue;
        found[count] = p.shape;
        count += 1;
    }
    return found[0..count];
}

/// How far `g`, placed by `xf`, has sunk into each shape it overlaps - or
/// come nearer to it than `filter.margin` - and the way out, a mesh
/// triangle by triangle. What a character steps out of what it has walked
/// into with.
pub fn penetrations(self: *World, g: Geometry, xf: Transform, filter: QueryFilter, found: []Penetration) []Penetration {
    var store: collide.Storage = .{};
    const mine = collide.solidOf(g, xf, &store) orelse return found[0..0];
    const Walk = struct {
        world: *World,
        mine: collide.Solid,
        box: Aabb,
        filter: QueryFilter,
        found: []Penetration,
        count: usize = 0,

        fn take(walk: *@This(), index: u32, m: Manifold) void {
            if (m.count == 0 or walk.count == walk.found.len) return;
            var deepest = m.points[0];
            for (m.pointSlice()) |p| if (p.separation < deepest.separation) {
                deepest = p;
            };
            if (deepest.separation >= walk.filter.margin) return;
            const entry = walk.world.shapeAt(index);
            walk.found[walk.count] = .{ .shape = walk.world.handleOfShape(index), .body = entry.body, .normal = m.normal, .depth = walk.filter.margin - deepest.separation, .point = deepest.point };
            walk.count += 1;
        }

        fn visit(walk: *@This(), index: u32) bool {
            const world = walk.world;
            const entry = world.shapeAt(index);
            if (!world.sees(entry, walk.filter)) return true;
            const shape_xf = world.shapeTransform(entry);
            if (entry.def.geometry == .mesh) {
                const mesh = entry.def.geometry.mesh;
                var triangles: TriangleList = .{ .gpa = world.gpa };
                defer triangles.list.deinit(world.gpa);
                mesh.query(boxInFrame(walk.box, shape_xf), &triangles, TriangleList.add);
                for (triangles.list.items) |t| {
                    const corners = mesh.triangle(t);
                    var store_tri: collide.Storage = .{};
                    const tri = collide.triangleSolid(&store_tri, shape_xf.apply(corners[0]), shape_xf.apply(corners[1]), shape_xf.apply(corners[2]));
                    walk.take(index, collide.collide(tri, walk.mine, walk.filter.margin));
                }
            } else {
                var other_store: collide.Storage = .{};
                const other = collide.solidOf(entry.def.geometry, shape_xf, &other_store).?;
                walk.take(index, collide.collide(other, walk.mine, walk.filter.margin));
            }
            return walk.count < walk.found.len;
        }
    };
    const box = g.aabb(xf);
    var walk: Walk = .{ .world = self, .mine = mine, .box = box, .filter = filter, .found = found };
    self.static_tree.query(box, &walk, Walk.visit);
    self.moving_tree.query(box, &walk, Walk.visit);
    return found[0..walk.count];
}

/// The shapes a point is inside.
pub fn overlapPoint(self: *World, point: Vec3, filter: QueryFilter, found: []ShapeId) []ShapeId {
    return self.overlapShape(.{ .sphere = .{ .radius = 0.001 } }, .at(point), filter, found);
}

/// The shapes whose boxes overlap `box`.
pub fn overlapAabb(self: *World, box: Aabb, filter: QueryFilter, found: []ShapeId) []ShapeId {
    const Walk = struct {
        world: *World,
        filter: QueryFilter,
        found: []ShapeId,
        count: usize = 0,
        box: Aabb,

        fn visit(walk: *@This(), index: u32) bool {
            const entry = walk.world.shapeAt(index);
            if (!walk.world.sees(entry, walk.filter)) return true;
            if (!entry.def.geometry.aabb(walk.world.shapeTransform(entry)).overlaps(walk.box)) return true;
            walk.found[walk.count] = walk.world.handleOfShape(index);
            walk.count += 1;
            return walk.count < walk.found.len;
        }
    };
    var walk: Walk = .{ .world = self, .filter = filter, .found = found, .box = box };
    self.static_tree.query(box, &walk, Walk.visit);
    if (walk.count < found.len) self.moving_tree.query(box, &walk, Walk.visit);
    return found[0..walk.count];
}

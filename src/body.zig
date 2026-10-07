// SPDX-License-Identifier: BSD-2-Clause

//! A rigid body: where it is, how it moves, and what it weighs.
//!
//! ```zig
//! const crate = try world.createBody(.{ .position = .init(0, 4, 0) });
//! _ = try world.addShape(crate, .box(.init(0.5, 0.5, 0.5)));
//! const b = world.body(crate).?;
//! b.applyImpulse(.init(0, 5, 0), b.center);
//! ```
//!
//! **Three kinds.** A *static* body never moves: the ground, the walls, the
//! level. A *kinematic* body moves the way it is told and nothing pushes
//! it: a lift, a door. A *dynamic* body falls, bounces, and is stopped by
//! the other two. Two bodies that cannot both move never collide.
//!
//! **It moves about its centre of mass, not its origin.** The origin is
//! where a game put it - a character's feet, a model's root - and the centre
//! of mass is where its shapes say it is. Velocity and turning act on the
//! centre, and the origin follows.
//!
//! **A dynamic body falls asleep** when it and everything it touches have
//! been still for a while, and a step passes it over until something wakes
//! it: a push, a new speed, a move, or something running into it.

const std = @import("std");
const id = @import("fluxion_id");

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Quat = geometry.Quat;
const Mat3 = geometry.Mat3;
const Transform = geometry.Transform;

const Body = @This();

pub const Id = id.Handle(Body);
pub const ShapeId = id.Handle(@import("World.zig").ShapeEntry);

pub const Type = enum(u8) {
    static,
    kinematic,
    dynamic,
};

/// What `World.createBody` takes.
pub const Def = struct {
    type: Type = .dynamic,
    position: Vec3 = .zero,
    rotation: Quat = .identity,
    linear_velocity: Vec3 = .zero,
    /// Radians a second about each axis.
    angular_velocity: Vec3 = .zero,
    /// How much of its speed it loses a second.
    linear_damping: f32 = 0,
    angular_damping: f32 = 0,
    /// Nought floats, minus one rises.
    gravity_scale: f32 = 1,
    /// It moves but never turns: a character standing up.
    lock_rotation: bool = false,
    allow_sleep: bool = true,
    /// The caller's: an entity, an index.
    user_data: u64 = 0,
};

type: Type,
/// Where its origin is, and how it is turned.
transform: Transform,
/// Its centre of mass, in the world and in its own frame.
center: Vec3,
local_center: Vec3 = .zero,
linear_velocity: Vec3,
angular_velocity: Vec3,
/// Pushing for the whole of the next step, then gone.
force: Vec3 = .zero,
torque: Vec3 = .zero,
mass: f32 = 0,
inv_mass: f32 = 0,
/// About its centre of mass, in its own frame.
inertia: Mat3 = .zero,
inv_inertia_local: Mat3 = .zero,
/// The same in the world, as it is turned now.
inv_inertia: Mat3 = .zero,
linear_damping: f32,
angular_damping: f32,
gravity_scale: f32,
lock_rotation: bool,
allow_sleep: bool,
awake: bool = true,
/// How long it has been nearly still.
sleep_time: f32 = 0,
first_shape: ShapeId = .none,
shape_count: u32 = 0,
user_data: u64,

pub fn init(def: Def) Body {
    const xf: Transform = .{ .position = def.position, .rotation = def.rotation.norm() };
    return .{
        .type = def.type,
        .transform = xf,
        .center = xf.position,
        .linear_velocity = if (def.type == .static) .zero else def.linear_velocity,
        .angular_velocity = if (def.type == .static) .zero else def.angular_velocity,
        .linear_damping = def.linear_damping,
        .angular_damping = def.angular_damping,
        .gravity_scale = def.gravity_scale,
        .lock_rotation = def.lock_rotation,
        .allow_sleep = def.allow_sleep,
        .user_data = def.user_data,
    };
}

pub fn position(self: *const Body) Vec3 {
    return self.transform.position;
}

pub fn rotation(self: *const Body) Quat {
    return self.transform.rotation;
}

/// How fast the point `p`, in the world, is moving.
pub fn velocityAt(self: *const Body, p: Vec3) Vec3 {
    return self.linear_velocity.add(self.angular_velocity.cross(p.sub(self.center)));
}

pub fn wake(self: *Body) void {
    if (self.type == .static) return;
    self.awake = true;
    self.sleep_time = 0;
}

pub fn setLinearVelocity(self: *Body, v: Vec3) void {
    if (self.type == .static) return;
    if (!v.eql(.zero)) self.wake();
    self.linear_velocity = v;
}

pub fn setAngularVelocity(self: *Body, w: Vec3) void {
    if (self.type == .static) return;
    if (!w.eql(.zero)) self.wake();
    self.angular_velocity = if (self.lock_rotation) .zero else w;
}

/// Push at `point` for the next step.
pub fn applyForce(self: *Body, f: Vec3, point: Vec3) void {
    if (self.type != .dynamic) return;
    self.wake();
    self.force = self.force.add(f);
    self.torque = self.torque.add(point.sub(self.center).cross(f));
}

pub fn applyForceToCenter(self: *Body, f: Vec3) void {
    if (self.type != .dynamic) return;
    self.wake();
    self.force = self.force.add(f);
}

pub fn applyTorque(self: *Body, t: Vec3) void {
    if (self.type != .dynamic) return;
    self.wake();
    self.torque = self.torque.add(t);
}

/// A kick at `point`: the speed changes at once.
pub fn applyImpulse(self: *Body, impulse: Vec3, point: Vec3) void {
    if (self.type != .dynamic) return;
    self.wake();
    self.linear_velocity = self.linear_velocity.mulAdd(impulse, self.inv_mass);
    self.angular_velocity = self.angular_velocity.add(geometry.mulVec(self.inv_inertia, point.sub(self.center).cross(impulse)));
}

pub fn applyImpulseToCenter(self: *Body, impulse: Vec3) void {
    if (self.type != .dynamic) return;
    self.wake();
    self.linear_velocity = self.linear_velocity.mulAdd(impulse, self.inv_mass);
}

pub fn applyAngularImpulse(self: *Body, impulse: Vec3) void {
    if (self.type != .dynamic) return;
    self.wake();
    self.angular_velocity = self.angular_velocity.add(geometry.mulVec(self.inv_inertia, impulse));
}

/// Its inverse inertia in the world, as it is turned now.
pub fn updateInertia(self: *Body) void {
    if (self.inv_mass == 0 or self.lock_rotation) {
        self.inv_inertia = .zero;
        return;
    }
    self.inv_inertia = geometry.rotateTensor(geometry.rotationOf(self.transform.rotation), self.inv_inertia_local);
}

/// Its origin where its centre of mass and its turn now say.
pub fn placeFromCenter(self: *Body) void {
    self.transform.position = self.center.sub(self.transform.rotation.rotate(self.local_center));
}

/// Whether it moves this step: awake and dynamic, or kinematic with
/// somewhere to go.
pub fn isActive(self: *const Body) bool {
    return switch (self.type) {
        .static => false,
        .kinematic => !self.linear_velocity.eql(.zero) or !self.angular_velocity.eql(.zero),
        .dynamic => self.awake,
    };
}

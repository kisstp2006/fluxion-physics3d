// SPDX-License-Identifier: BSD-2-Clause

//! Joints: two bodies held to each other at a point, turning only as far as
//! the joint lets them.
//!
//! ```zig
//! // A shoulder: the arm pinned to the chest, swinging in a cone of half a
//! // right angle round where it points, and twisting a little each way.
//! _ = try world.createJoint(.{
//!     .body_a = chest,
//!     .body_b = arm,
//!     .frame_a = .{ .position = .init(0.2, 0.4, 0) },
//!     .frame_b = .{},
//!     .kind = .ball,
//!     .swing = std.math.pi / 4.0,
//!     .twist = .{ .lower = -0.3, .upper = 0.3 },
//! });
//!
//! // A knee: turning about one axis only, from straight to bent.
//! _ = try world.createJoint(.{ .body_a = thigh, .body_b = shin, .frame_a = knee_on_thigh, .frame_b = knee_on_shin, .kind = .hinge, .twist = .{ .lower = 0, .upper = 2 } });
//! ```
//!
//! | Joint | What it holds | What it limits |
//! | --- | --- | --- |
//! | `ball` | the two frames' origins together | how far the frames' `x` axes swing apart - a cone - and how far one twists about it |
//! | `hinge` | the origins together, and the frames' `z` axes in line | how far one turns about `z` |
//!
//! **Each body has a frame of the joint on it**, given in the body's own
//! space: where the joint is, and how it is turned there. The joint is at
//! rest where the two frames are one. A ragdoll's joints are made from its
//! bones as they were drawn, not as they stand when it falls: the frames
//! are what says where its rest is, whatever pose the bodies are in when
//! the joint is made. `framesAt` gives the two for a joint placed in the
//! world, with the bodies' pose now its rest.
//!
//! **A ball's limits, a cone and a twist**, as a shoulder has them: `swing`
//! is how far the `x` axis of `frame_b` may lean from that of `frame_a`,
//! any way round; `twist` how far `frame_b` may turn about its `x` axis,
//! the swing taken out, from `lower` to `upper`. A hinge's `twist` is its
//! angle about `z`. A limit is like a contact: it pushes only when it is
//! reached, and a joint coming up to one fast stops at it within the
//! substep.
//!
//! **Solved with the contacts, and the same way**: sequential impulses,
//! warm-started from what the same joint pushed the step before, measured
//! again from where the bodies are at every pass, drift taken back as a
//! stiff, heavily damped spring would take it (`Softness`) and then
//! relaxed out of the velocities. The joints go first in every pass.
//!
//! **Two bodies a joint holds do not collide with each other** unless
//! `collide` says so: a joint usually holds two shapes that overlap on
//! purpose - an arm in its shoulder - and the contact would fight it.

const std = @import("std");
const testing = std.testing;
const id = @import("fluxion_id");

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Quat = geometry.Quat;
const Mat3 = geometry.Mat3;
const Transform = geometry.Transform;
const mulVec = geometry.mulVec;
const Body = @import("body.zig");
const Softness = @import("Softness.zig");

pub const Id = id.Handle(Joint);

pub const Kind = enum(u8) { ball, hinge };

/// A range of turn, in radians.
pub const Range = struct {
    lower: f32 = -std.math.pi,
    upper: f32 = std.math.pi,
};

pub const Def = struct {
    body_a: Body.Id,
    body_b: Body.Id,
    /// The joint on each body, in its own space.
    frame_a: Transform = .identity,
    frame_b: Transform = .identity,
    kind: Kind = .ball,
    /// A ball's cone: how far its `x` axes may lean apart. Null leaves the
    /// swing free.
    swing: ?f32 = null,
    /// A ball's turn about `x`, or a hinge's about `z`. Null leaves it free.
    twist: ?Range = null,
    /// Whether the two bodies still collide with each other.
    collide: bool = false,
};

/// The two frames of a joint placed at `at`, in the world, on two bodies
/// placed at `a` and `b`: the pose they are in now its rest.
pub fn framesAt(a: Transform, b: Transform, at: Transform) [2]Transform {
    return .{ a.invMul(at), b.invMul(at) };
}

pub const Joint = struct {
    body_a: Body.Id,
    body_b: Body.Id,
    frame_a: Transform,
    frame_b: Transform,
    kind: Kind,
    swing: ?f32,
    twist: ?Range,
    collide: bool,
    /// What it pushed, kept from one pass to the next and from the step
    /// before: the point's, the hinge's two that keep its axes in line,
    /// the swing's, and the twist's at each end.
    point_impulse: Vec3 = .zero,
    align_impulse: [2]f32 = .{ 0, 0 },
    swing_impulse: f32 = 0,
    twist_impulse: [2]f32 = .{ 0, 0 },

    pub fn init(def: Def) Joint {
        return .{
            .body_a = def.body_a,
            .body_b = def.body_b,
            .frame_a = .{ .position = def.frame_a.position, .rotation = def.frame_a.rotation.norm() },
            .frame_b = .{ .position = def.frame_b.position, .rotation = def.frame_b.rotation.norm() },
            .kind = def.kind,
            .swing = def.swing,
            .twist = def.twist,
            .collide = def.collide,
        };
    }

    /// Its two frames in the world, as the bodies are now.
    pub fn worldFrames(self: *const Joint, a: *const Body, b: *const Body) [2]Transform {
        return .{ a.transform.mul(self.frame_a), b.transform.mul(self.frame_b) };
    }

    /// How far `frame_b` is turned about the twist's axis from `frame_a`,
    /// the swing taken out: in (-pi, pi].
    pub fn twistAngle(self: *const Joint, a: *const Body, b: *const Body) f32 {
        const frames = self.worldFrames(a, b);
        return twistOf(frames[0].rotation.inverse().mul(frames[1].rotation), self.kind);
    }

    /// How far the `x` axes lean apart.
    pub fn swingAngle(self: *const Joint, a: *const Body, b: *const Body) f32 {
        const frames = self.worldFrames(a, b);
        const xa = frames[0].rotation.rotate(.unit_x);
        const xb = frames[1].rotation.rotate(.unit_x);
        return std.math.acos(std.math.clamp(xa.dot(xb), -1, 1));
    }
};

/// The turn of `relative` about a ball's `x` or a hinge's `z`, the rest of
/// it taken out.
pub fn twistOf(relative: Quat, kind: Kind) f32 {
    const along = if (kind == .ball) relative.x else relative.z;
    var angle = 2 * std.math.atan2(along, relative.w);
    if (angle > std.math.pi) angle -= 2 * std.math.pi;
    if (angle <= -std.math.pi) angle += 2 * std.math.pi;
    return angle;
}

pub const Pass = enum { solve, relax };

/// What every pass is told, worked out once a step.
pub const Step = struct {
    inv_h: f32,
    softness: Softness,
    /// The fastest a limit past its end is turned back, in radians a
    /// second: a ragdoll that starts in a pose past its limits - an arm
    /// raised to aim - swings back into them, and does not snap.
    max_turn: f32,
};

/// `r` as the matrix that crosses with it: `skew(r) v = r x v`.
fn skew(r: Vec3) Mat3 {
    return .fromCols(.init(0, r.z, -r.y), .init(-r.z, 0, r.x), .init(r.y, -r.x, 0));
}

/// The two bodies' velocities, worked on locally and written back.
const Velocities = struct {
    va: Vec3,
    wa: Vec3,
    vb: Vec3,
    wb: Vec3,

    fn of(a: *const Body, b: *const Body) Velocities {
        return .{ .va = a.linear_velocity, .wa = a.angular_velocity, .vb = b.linear_velocity, .wb = b.angular_velocity };
    }

    fn pushPoint(v: *Velocities, a: *const Body, b: *const Body, impulse: Vec3, ra: Vec3, rb: Vec3) void {
        v.va = v.va.mulAdd(impulse, -a.inv_mass);
        v.wa = v.wa.sub(mulVec(a.inv_inertia, ra.cross(impulse)));
        v.vb = v.vb.mulAdd(impulse, b.inv_mass);
        v.wb = v.wb.add(mulVec(b.inv_inertia, rb.cross(impulse)));
    }

    fn pushTurn(v: *Velocities, a: *const Body, b: *const Body, impulse: Vec3) void {
        v.wa = v.wa.sub(mulVec(a.inv_inertia, impulse));
        v.wb = v.wb.add(mulVec(b.inv_inertia, impulse));
    }

    fn store(v: Velocities, a: *Body, b: *Body) void {
        if (a.type == .dynamic) {
            a.linear_velocity = v.va;
            a.angular_velocity = if (a.lock_rotation) .zero else v.wa;
        }
        if (b.type == .dynamic) {
            b.linear_velocity = v.vb;
            b.angular_velocity = if (b.lock_rotation) .zero else v.wb;
        }
    }
};

/// What one unit of turning impulse about `n` changes the two bodies'
/// relative spin about it by.
fn turnMass(a: *const Body, b: *const Body, n: Vec3) f32 {
    const k = n.dot(mulVec(a.inv_inertia, n)) + n.dot(mulVec(b.inv_inertia, n));
    return if (k > 0) 1 / k else 0;
}

/// Push again what the joint pushed before.
pub fn warmStart(j: *const Joint, a: *Body, b: *Body) void {
    var v: Velocities = .of(a, b);
    const frames = j.worldFrames(a, b);
    const ra = frames[0].position.sub(a.center);
    const rb = frames[1].position.sub(b.center);
    v.pushPoint(a, b, j.point_impulse, ra, rb);
    const axes = turnAxes(j, frames);
    var turn = axes.swing.scale(-j.swing_impulse).add(axes.twist.scale(j.twist_impulse[0] - j.twist_impulse[1]));
    if (j.kind == .hinge) turn = turn.add(axes.perpendicular[0].scale(j.align_impulse[0])).add(axes.perpendicular[1].scale(j.align_impulse[1]));
    v.pushTurn(a, b, turn);
    v.store(a, b);
}

const Axes = struct {
    /// The twist's: a ball's `x`, between the two bodies' own; a hinge's `z`.
    twist: Vec3,
    /// Which way leaning further apart turns, for a ball's cone.
    swing: Vec3,
    /// A hinge's: at right angles to its axis.
    perpendicular: [2]Vec3,
};

fn turnAxes(j: *const Joint, frames: [2]Transform) Axes {
    switch (j.kind) {
        .ball => {
            const xa = frames[0].rotation.rotate(.unit_x);
            const xb = frames[1].rotation.rotate(.unit_x);
            const lean = xa.cross(xb);
            const lean_len = lean.len();
            const between = xa.add(xb);
            return .{
                .twist = if (between.lenSq() > 1e-8) between.norm() else xb,
                .swing = if (lean_len > 1e-6) lean.scale(1 / lean_len) else .zero,
                .perpendicular = .{ .zero, .zero },
            };
        },
        .hinge => {
            const za = frames[0].rotation.rotate(.unit_z);
            return .{ .twist = za, .swing = .zero, .perpendicular = geometry.basis(za) };
        },
    }
}

/// A limit, solved as a contact's normal is: `gap` is how far it is from
/// being reached, and it pushes along `n` only to keep it from closing
/// past nought.
fn limit(v: *Velocities, a: *const Body, b: *const Body, n: Vec3, gap: f32, total: *f32, step: Step, comptime pass: Pass) void {
    const mass = turnMass(a, b, n);
    if (mass == 0) return;
    var bias: f32 = 0;
    var mass_scale: f32 = 1;
    var impulse_scale: f32 = 0;
    if (gap > 0) {
        bias = gap * step.inv_h;
    } else if (pass == .solve) {
        bias = @max(step.softness.bias_rate * gap, -step.max_turn);
        mass_scale = step.softness.mass_scale;
        impulse_scale = step.softness.impulse_scale;
    }
    const speed = v.wb.sub(v.wa).dot(n);
    const lambda = -mass * mass_scale * (speed + bias) - impulse_scale * total.*;
    const kept = @max(total.* + lambda, 0);
    const applied = kept - total.*;
    total.* = kept;
    v.pushTurn(a, b, n.scale(applied));
}

/// One pass over the joint: what keeps it together, then its limits.
pub fn solve(j: *Joint, a: *Body, b: *Body, step: Step, comptime pass: Pass) void {
    var v: Velocities = .of(a, b);
    const frames = j.worldFrames(a, b);
    const soft: Softness = if (pass == .solve) step.softness else .rigid;

    // A hinge's axes kept in line.
    if (j.kind == .hinge) {
        const za = frames[0].rotation.rotate(.unit_z);
        const zb = frames[1].rotation.rotate(.unit_z);
        const err = za.cross(zb);
        for (geometry.basis(za), 0..) |p, k| {
            const mass = turnMass(a, b, p);
            if (mass == 0) continue;
            const speed = v.wb.sub(v.wa).dot(p);
            const lambda = soft.impulse(mass, speed, err.dot(p), j.align_impulse[k]);
            j.align_impulse[k] += lambda;
            v.pushTurn(a, b, p.scale(lambda));
        }
    }

    // The limits.
    const axes = turnAxes(j, frames);
    if (j.kind == .ball) if (j.swing) |cone| {
        if (!axes.swing.eql(.zero)) {
            const xa = frames[0].rotation.rotate(.unit_x);
            const xb = frames[1].rotation.rotate(.unit_x);
            const lean = std.math.acos(std.math.clamp(xa.dot(xb), -1, 1));
            limit(&v, a, b, axes.swing.neg(), cone - lean, &j.swing_impulse, step, pass);
        }
    };
    if (j.twist) |range| {
        const turned = twistOf(frames[0].rotation.inverse().mul(frames[1].rotation), j.kind);
        limit(&v, a, b, axes.twist, turned - range.lower, &j.twist_impulse[0], step, pass);
        limit(&v, a, b, axes.twist.neg(), range.upper - turned, &j.twist_impulse[1], step, pass);
    }

    // The two frames' origins together, all three ways at once.
    const ra = frames[0].position.sub(a.center);
    const rb = frames[1].position.sub(b.center);
    const sa = skew(ra);
    const sb = skew(rb);
    // K = (ma + mb) I - [ra] Ia [ra] - [rb] Ib [rb]
    var k = geometry.diagonal(.splat(a.inv_mass + b.inv_mass));
    k = geometry.addMat(k, geometry.scaleMat(geometry.mulMat(geometry.mulMat(sa, a.inv_inertia), sa), -1));
    k = geometry.addMat(k, geometry.scaleMat(geometry.mulMat(geometry.mulMat(sb, b.inv_inertia), sb), -1));
    const mass = geometry.inverseOrZero(k);
    const separation = frames[1].position.sub(frames[0].position);
    const speed = v.vb.add(v.wb.cross(rb)).sub(v.va).sub(v.wa.cross(ra));
    const wanted = mulVec(mass, speed.add(separation.scale(soft.bias_rate))).scale(-soft.mass_scale).sub(j.point_impulse.scale(soft.impulse_scale));
    j.point_impulse = j.point_impulse.add(wanted);
    v.pushPoint(a, b, wanted, ra, rb);
    v.store(a, b);
}

test "the twist of a turn about x, or about z for a hinge, with the rest taken out" {
    const quarter = std.math.pi / 2.0;
    try testing.expectApproxEqAbs(quarter, twistOf(.fromAxisAngle(.unit_x, quarter), .ball), 1e-5);
    try testing.expectApproxEqAbs(@as(f32, 0), twistOf(.fromAxisAngle(.unit_y, quarter), .ball), 1e-5);
    try testing.expectApproxEqAbs(-quarter, twistOf(.fromAxisAngle(.unit_z, -quarter), .hinge), 1e-5);
    // A swing then a twist: the twist alone is read.
    const both = Quat.fromAxisAngle(.unit_y, 0.4).mul(Quat.fromAxisAngle(.unit_x, 0.3));
    try testing.expectApproxEqAbs(@as(f32, 0.3), twistOf(both, .ball), 1e-5);
}

test "the frames of a joint placed in the world are where it is, seen from each body" {
    const a: Transform = .{ .position = .init(1, 0, 0), .rotation = .fromAxisAngle(.unit_y, 0.5) };
    const b: Transform = .{ .position = .init(0, 2, 0) };
    const at: Transform = .{ .position = .init(0.5, 1, 0), .rotation = .fromAxisAngle(.unit_z, 0.25) };
    const frames = framesAt(a, b, at);
    const back_a = a.mul(frames[0]);
    const back_b = b.mul(frames[1]);
    try testing.expect(back_a.position.sub(at.position).len() < 1e-5);
    try testing.expect(back_b.position.sub(at.position).len() < 1e-5);
    try testing.expect(back_a.rotation.approxEql(at.rotation));
}

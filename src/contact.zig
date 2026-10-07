// SPDX-License-Identifier: BSD-2-Clause

//! A contact as the solver sees it: two bodies, a normal, and for each
//! point the numbers that turn a relative velocity into an impulse.
//!
//! **Sequential impulses, in substeps.** Each contact is solved on its own,
//! each pass correcting what the others undid, and a step is cut into a few
//! substeps that each integrate the velocities, solve once, move the bodies
//! and relax once. The manifold is found once a step, but how deep each
//! point has sunk is worked out again every substep from where the bodies
//! are now: each point is kept in both bodies' own frames, and the gap
//! between the two copies is the depth. This is Erin Catto's soft step.
//!
//! **Sinking in is taken back as a stiff, damped spring would take it** -
//! see `Softness` - never faster than a limit, and a relaxing pass with no
//! push at all leaves the velocities what they would have been had nothing
//! needed pushing apart. A point still apart may close the gap by the end of
//! the substep and no more: a speculative contact, which stops a fast body
//! at a surface before it can pass through it.
//!
//! **Friction pulls two ways along the surface**, no harder in all than the
//! friction times how hard the surfaces press together: a cone, not the
//! square a pull along each way alone would make.
//!
//! **Restitution comes last**, from the speed each point came in at, and
//! only where it pushed.

const std = @import("std");

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Mat3 = geometry.Mat3;
const mulVec = geometry.mulVec;
const Body = @import("body.zig");
const collide = @import("collide.zig");
const Manifold = collide.Manifold;
const Softness = @import("Softness.zig");

/// What a contact is remembered by from one step to the next: its two
/// shapes, as their handles' bits, and which triangle for a mesh.
pub const PairKey = struct {
    a: u64,
    b: u64,
    part: u32 = 0,
};

/// What a point ended a step with, kept for the next one.
pub const Impulse = struct {
    id: u32 = 0,
    normal: f32 = 0,
    tangent: [2]f32 = .{ 0, 0 },
};

/// What every pass is told, worked out once a step.
pub const Step = struct {
    /// One over the substep.
    inv_h: f32,
    /// How deep a point may sit before it is pushed.
    slop: f32,
    /// The fastest a push may be.
    max_push: f32,
};

pub const PointConstraint = struct {
    /// From each body's centre of mass to the point, at the start of the
    /// step.
    ra: Vec3,
    rb: Vec3,
    /// The point in each body's own frame.
    local_a: Vec3,
    local_b: Vec3,
    separation: f32,
    normal_mass: f32,
    tangent_mass: [2]f32,
    /// How fast the two were coming together here before anything was
    /// solved. Negative when approaching.
    approach: f32,
    normal_impulse: f32,
    tangent_impulse: [2]f32,
    max_normal_impulse: f32,
    id: u32,
};

pub const Constraint = struct {
    key: PairKey,
    body_a: u32,
    body_b: u32,
    normal: Vec3,
    tangents: [2]Vec3,
    friction: f32,
    restitution: f32,
    inv_mass_a: f32,
    inv_mass_b: f32,
    inv_inertia_a: Mat3,
    inv_inertia_b: Mat3,
    softness: Softness,
    points: [collide.max_points]PointConstraint,
    count: u32,

    pub fn pointSlice(self: *const Constraint) []const PointConstraint {
        return self.points[0..self.count];
    }
};

/// What one unit of impulse along `n` at `r` from the centre changes the
/// speed there by, for one body.
inline fn angularPart(inv_inertia: Mat3, r: Vec3, n: Vec3) f32 {
    const rn = r.cross(n);
    return rn.dot(mulVec(inv_inertia, rn));
}

/// The constraint for one manifold. `warm` is what the same pair ended the
/// last step with, or null.
pub fn prepare(
    key: PairKey,
    manifold: *const Manifold,
    a: *const Body,
    b: *const Body,
    index_a: u32,
    index_b: u32,
    friction: f32,
    restitution: f32,
    warm: ?[collide.max_points]Impulse,
    softness: Softness,
    static_softness: Softness,
) Constraint {
    const n = manifold.normal;
    const t = geometry.basis(n);
    var c: Constraint = .{
        .key = key,
        .body_a = index_a,
        .body_b = index_b,
        .normal = n,
        .tangents = t,
        .friction = friction,
        .restitution = restitution,
        .inv_mass_a = a.inv_mass,
        .inv_mass_b = b.inv_mass,
        .inv_inertia_a = a.inv_inertia,
        .inv_inertia_b = b.inv_inertia,
        .softness = if (a.inv_mass == 0 or b.inv_mass == 0) static_softness else softness,
        .points = undefined,
        .count = manifold.count,
    };
    const masses = c.inv_mass_a + c.inv_mass_b;
    for (manifold.pointSlice(), 0..) |mp, i| {
        const ra = mp.point.sub(a.center);
        const rb = mp.point.sub(b.center);
        const k_normal = masses + angularPart(c.inv_inertia_a, ra, n) + angularPart(c.inv_inertia_b, rb, n);
        var tangent_mass: [2]f32 = undefined;
        for (0..2) |k| {
            const k_tangent = masses + angularPart(c.inv_inertia_a, ra, t[k]) + angularPart(c.inv_inertia_b, rb, t[k]);
            tangent_mass[k] = if (k_tangent > 0) 1 / k_tangent else 0;
        }
        const dv = b.velocityAt(mp.point).sub(a.velocityAt(mp.point));
        var point: PointConstraint = .{
            .ra = ra,
            .rb = rb,
            .local_a = a.transform.unapply(mp.point),
            .local_b = b.transform.unapply(mp.point),
            .separation = mp.separation,
            .normal_mass = if (k_normal > 0) 1 / k_normal else 0,
            .tangent_mass = tangent_mass,
            .approach = dv.dot(n),
            .normal_impulse = 0,
            .tangent_impulse = .{ 0, 0 },
            .max_normal_impulse = 0,
            .id = mp.id,
        };
        if (warm) |last| {
            for (last) |w| {
                if (w.id == mp.id and (w.normal != 0 or w.tangent[0] != 0 or w.tangent[1] != 0)) {
                    point.normal_impulse = w.normal;
                    point.tangent_impulse = w.tangent;
                    break;
                }
            }
        }
        c.points[i] = point;
    }
    return c;
}

pub const Pass = enum { solve, relax };

/// The velocities of the two bodies, copied out so a pass works on locals,
/// and written back after.
const Velocities = struct {
    va: Vec3,
    wa: Vec3,
    vb: Vec3,
    wb: Vec3,

    fn of(a: *const Body, b: *const Body) Velocities {
        return .{ .va = a.linear_velocity, .wa = a.angular_velocity, .vb = b.linear_velocity, .wb = b.angular_velocity };
    }

    inline fn along(v: *const Velocities, n: Vec3, ra: Vec3, rb: Vec3) f32 {
        return v.vb.add(v.wb.cross(rb)).sub(v.va).sub(v.wa.cross(ra)).dot(n);
    }

    inline fn push(v: *Velocities, c: *const Constraint, impulse: Vec3, ra: Vec3, rb: Vec3) void {
        v.va = v.va.mulAdd(impulse, -c.inv_mass_a);
        v.wa = v.wa.sub(mulVec(c.inv_inertia_a, ra.cross(impulse)));
        v.vb = v.vb.mulAdd(impulse, c.inv_mass_b);
        v.wb = v.wb.add(mulVec(c.inv_inertia_b, rb.cross(impulse)));
    }

    fn store(v: Velocities, c: *const Constraint, a: *Body, b: *Body) void {
        if (c.inv_mass_a != 0) {
            a.linear_velocity = v.va;
            a.angular_velocity = v.wa;
        }
        if (c.inv_mass_b != 0) {
            b.linear_velocity = v.vb;
            b.angular_velocity = v.wb;
        }
    }
};

pub fn warmStart(c: *const Constraint, a: *Body, b: *Body) void {
    var v: Velocities = .of(a, b);
    for (c.pointSlice()) |p| {
        const impulse = c.normal.scale(p.normal_impulse).add(c.tangents[0].scale(p.tangent_impulse[0])).add(c.tangents[1].scale(p.tangent_impulse[1]));
        v.push(c, impulse, p.ra, p.rb);
    }
    v.store(c, a, b);
}

/// How deep a point is now: its depth when the manifold was made, and how
/// far its two copies - one on each body - have moved apart along the
/// normal since. Negative is sunk in.
inline fn separationNow(c: *const Constraint, p: *const PointConstraint, a: *const Body, b: *const Body) f32 {
    const pa = a.transform.apply(p.local_a);
    const pb = b.transform.apply(p.local_b);
    return p.separation + pb.sub(pa).dot(c.normal);
}

/// One pass: the normal of each point, then its friction - which can be no
/// more than how hard the normal pressed.
pub fn solve(c: *Constraint, a: *Body, b: *Body, step: Step, comptime pass: Pass) void {
    var v: Velocities = .of(a, b);
    for (c.points[0..c.count]) |*p| {
        const s = separationNow(c, p, a, b);
        var bias: f32 = 0;
        var mass_scale: f32 = 1;
        var impulse_scale: f32 = 0;
        if (s > 0) {
            // Apart: it may close the gap by the end of the substep.
            bias = s * step.inv_h;
        } else if (pass == .solve) {
            bias = @max(c.softness.bias_rate * @min(0, s + step.slop), -step.max_push);
            mass_scale = c.softness.mass_scale;
            impulse_scale = c.softness.impulse_scale;
        }
        const vn = v.along(c.normal, p.ra, p.rb);
        var lambda = -p.normal_mass * mass_scale * (vn + bias) - impulse_scale * p.normal_impulse;
        // It pushes and never pulls: the total is clamped, not this pass's.
        const total = @max(p.normal_impulse + lambda, 0);
        lambda = total - p.normal_impulse;
        p.normal_impulse = total;
        p.max_normal_impulse = @max(p.max_normal_impulse, lambda);
        v.push(c, c.normal.scale(lambda), p.ra, p.rb);
    }
    for (c.points[0..c.count]) |*p| {
        var wanted: [2]f32 = undefined;
        for (0..2) |k| {
            const vt = v.along(c.tangents[k], p.ra, p.rb);
            wanted[k] = p.tangent_impulse[k] - p.tangent_mass[k] * vt;
        }
        // Inside the cone: no more than the friction times the push.
        const most = c.friction * p.normal_impulse;
        const length = @sqrt(wanted[0] * wanted[0] + wanted[1] * wanted[1]);
        if (length > most) {
            const shrink = if (length > 0) most / length else 0;
            wanted[0] *= shrink;
            wanted[1] *= shrink;
        }
        const d0 = wanted[0] - p.tangent_impulse[0];
        const d1 = wanted[1] - p.tangent_impulse[1];
        p.tangent_impulse = wanted;
        v.push(c, c.tangents[0].scale(d0).add(c.tangents[1].scale(d1)), p.ra, p.rb);
    }
    v.store(c, a, b);
}

/// A bounce from the speed each point came in at, where it came in fast
/// enough and pushed.
pub fn restitute(c: *Constraint, a: *Body, b: *Body, threshold: f32) void {
    if (c.restitution == 0) return;
    var v: Velocities = .of(a, b);
    for (c.points[0..c.count]) |*p| {
        if (p.approach > -threshold or p.max_normal_impulse == 0) continue;
        const vn = v.along(c.normal, p.ra, p.rb);
        var lambda = -p.normal_mass * (vn + c.restitution * p.approach);
        const total = @max(p.normal_impulse + lambda, 0);
        lambda = total - p.normal_impulse;
        p.normal_impulse = total;
        v.push(c, c.normal.scale(lambda), p.ra, p.rb);
    }
    v.store(c, a, b);
}

/// What each point ended the step with, for the next one.
pub fn impulsesOf(c: *const Constraint) [collide.max_points]Impulse {
    var out: [collide.max_points]Impulse = @splat(.{});
    for (c.pointSlice(), 0..) |p, i| out[i] = .{ .id = p.id, .normal = p.normal_impulse, .tangent = p.tangent_impulse };
    return out;
}

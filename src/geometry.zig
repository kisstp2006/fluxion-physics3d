// SPDX-License-Identifier: BSD-2-Clause

//! Places, turns and boxes in three dimensions, and the few matrix sums the
//! rest needs.
//!
//! **The axes are the engine's**: `+y` up, a right-handed frame, so `+z`
//! comes out of the screen towards whoever looks down `-z`. Gravity pulls
//! towards `-y`. A turn is a unit quaternion, kept unit by renormalising
//! after each step.

const std = @import("std");
const testing = std.testing;
const math = @import("fluxion_math");

pub const Vec3 = math.Vec3;
pub const Quat = math.Quat;
pub const Mat3 = math.Mat3;

/// Where something is and how it is turned: a point moves by the turn and
/// then the place. No scale - a shape's size is its own.
pub const Transform = struct {
    position: Vec3 = .zero,
    rotation: Quat = .identity,

    pub const identity: Transform = .{};

    pub fn at(position: Vec3) Transform {
        return .{ .position = position };
    }

    /// A point of this frame, in the frame outside it.
    pub inline fn apply(self: Transform, p: Vec3) Vec3 {
        return self.rotation.rotate(p).add(self.position);
    }

    /// A point of the frame outside, in this one.
    pub inline fn unapply(self: Transform, p: Vec3) Vec3 {
        return self.rotation.unrotate(p.sub(self.position));
    }

    /// A direction of this frame, outside it.
    pub inline fn turn(self: Transform, v: Vec3) Vec3 {
        return self.rotation.rotate(v);
    }

    pub inline fn unturn(self: Transform, v: Vec3) Vec3 {
        return self.rotation.unrotate(v);
    }

    /// `b` within `a`: what applying `b` and then `a` does.
    pub fn mul(a: Transform, b: Transform) Transform {
        return .{ .position = a.apply(b.position), .rotation = a.rotation.mul(b.rotation) };
    }

    /// `b` as seen from `a`: what `a.inverse().mul(b)` is, in one go.
    pub fn invMul(a: Transform, b: Transform) Transform {
        return .{
            .position = a.rotation.unrotate(b.position.sub(a.position)),
            .rotation = a.rotation.inverse().mul(b.rotation),
        };
    }

    pub fn inverse(self: Transform) Transform {
        const back = self.rotation.inverse();
        return .{ .position = back.rotate(self.position).neg(), .rotation = back };
    }
};

/// A box along the axes.
pub const Aabb = struct {
    min: Vec3,
    max: Vec3,

    /// Nothing: joined with anything, that thing.
    pub const empty: Aabb = .{ .min = .splat(std.math.inf(f32)), .max = .splat(-std.math.inf(f32)) };

    pub fn ofPoints(points: []const Vec3) Aabb {
        var box = empty;
        for (points) |p| box = box.add(p);
        return box;
    }

    /// The box that holds this one and `p`.
    pub inline fn add(self: Aabb, p: Vec3) Aabb {
        return .{ .min = self.min.min(p), .max = self.max.max(p) };
    }

    pub inline fn join(a: Aabb, b: Aabb) Aabb {
        return .{ .min = a.min.min(b.min), .max = a.max.max(b.max) };
    }

    /// Bigger by `margin` on every side.
    pub inline fn grow(self: Aabb, margin: f32) Aabb {
        return .{ .min = self.min.sub(.splat(margin)), .max = self.max.add(.splat(margin)) };
    }

    /// Stretched along `motion` to hold where it would get to.
    pub fn sweep(self: Aabb, motion: Vec3) Aabb {
        return .{ .min = self.min.add(motion.min(.zero)), .max = self.max.add(motion.max(.zero)) };
    }

    pub inline fn overlaps(a: Aabb, b: Aabb) bool {
        return a.min.x <= b.max.x and a.max.x >= b.min.x and
            a.min.y <= b.max.y and a.max.y >= b.min.y and
            a.min.z <= b.max.z and a.max.z >= b.min.z;
    }

    pub inline fn contains(self: Aabb, inner: Aabb) bool {
        return self.min.x <= inner.min.x and self.min.y <= inner.min.y and self.min.z <= inner.min.z and
            self.max.x >= inner.max.x and self.max.y >= inner.max.y and self.max.z >= inner.max.z;
    }

    pub inline fn containsPoint(self: Aabb, p: Vec3) bool {
        return p.x >= self.min.x and p.y >= self.min.y and p.z >= self.min.z and
            p.x <= self.max.x and p.y <= self.max.y and p.z <= self.max.z;
    }

    pub inline fn center(self: Aabb) Vec3 {
        return self.min.add(self.max).scale(0.5);
    }

    /// Half its size.
    pub inline fn extent(self: Aabb) Vec3 {
        return self.max.sub(self.min).scale(0.5);
    }

    /// What a tree of boxes keeps small: the chance a random ray or box
    /// meets one is in proportion to its surface.
    pub inline fn surfaceArea(self: Aabb) f32 {
        const d = self.max.sub(self.min);
        return 2 * (d.x * d.y + d.y * d.z + d.z * d.x);
    }

    /// Where a ray from `origin` along `translation` first enters the box,
    /// as a fraction of `translation` no more than `max_fraction`, or null
    /// for a miss. A ray that starts inside enters at nought.
    pub fn rayCast(self: Aabb, origin: Vec3, translation: Vec3, max_fraction: f32) ?f32 {
        var low: f32 = 0;
        var high: f32 = max_fraction;
        const o = [3]f32{ origin.x, origin.y, origin.z };
        const d = [3]f32{ translation.x, translation.y, translation.z };
        const lo = [3]f32{ self.min.x, self.min.y, self.min.z };
        const hi = [3]f32{ self.max.x, self.max.y, self.max.z };
        for (0..3) |axis| {
            if (@abs(d[axis]) < 1e-12) {
                if (o[axis] < lo[axis] or o[axis] > hi[axis]) return null;
                continue;
            }
            const inv = 1 / d[axis];
            var t1 = (lo[axis] - o[axis]) * inv;
            var t2 = (hi[axis] - o[axis]) * inv;
            if (t1 > t2) std.mem.swap(f32, &t1, &t2);
            low = @max(low, t1);
            high = @min(high, t2);
            if (low > high) return null;
        }
        return low;
    }
};

// -------------------------------------------------------------------------
// Matrices
// -------------------------------------------------------------------------

/// `m` times the column `v`.
pub inline fn mulVec(m: Mat3, v: Vec3) Vec3 {
    return m.cols[0].scale(v.x).add(m.cols[1].scale(v.y)).add(m.cols[2].scale(v.z));
}

pub fn transpose(m: Mat3) Mat3 {
    return .fromRows(m.cols[0], m.cols[1], m.cols[2]);
}

pub fn mulMat(a: Mat3, b: Mat3) Mat3 {
    return .fromCols(mulVec(a, b.cols[0]), mulVec(a, b.cols[1]), mulVec(a, b.cols[2]));
}

pub fn addMat(a: Mat3, b: Mat3) Mat3 {
    return .fromCols(a.cols[0].add(b.cols[0]), a.cols[1].add(b.cols[1]), a.cols[2].add(b.cols[2]));
}

pub fn scaleMat(m: Mat3, s: f32) Mat3 {
    return .fromCols(m.cols[0].scale(s), m.cols[1].scale(s), m.cols[2].scale(s));
}

pub fn diagonal(d: Vec3) Mat3 {
    return .fromScale(d);
}

/// `a b^T`.
pub fn outer(a: Vec3, b: Vec3) Mat3 {
    return .fromCols(a.scale(b.x), a.scale(b.y), a.scale(b.z));
}

/// The rotation `q` names, as a matrix.
pub fn rotationOf(q: Quat) Mat3 {
    return .fromQuat(q);
}

/// `r m r^T`: a tensor of a frame turned by `r`, in the frame outside it.
pub fn rotateTensor(r: Mat3, m: Mat3) Mat3 {
    return mulMat(mulMat(r, m), transpose(r));
}

/// The inverse of `m`, or zero where it has none: a body that cannot turn
/// about some axis.
pub fn inverseOrZero(m: Mat3) Mat3 {
    const det = m.det();
    if (@abs(det) < 1e-20) return .zero;
    return m.inverse() orelse .zero;
}

/// Two unit vectors at right angles to `n` and to each other, making a
/// right-handed frame with it: the two ways friction pulls.
pub fn basis(n: Vec3) [2]Vec3 {
    // From whichever axis `n` is least along, so the cross is never small.
    const t1 = if (@abs(n.x) >= 0.57735) Vec3.init(n.y, -n.x, 0).norm() else Vec3.init(0, n.z, -n.y).norm();
    return .{ t1, n.cross(t1) };
}

/// The closest points of two segments `p1 q1` and `p2 q2`, as fractions
/// along each: Ericson's `ClosestPtSegmentSegment`.
pub fn closestSegmentSegment(p1: Vec3, q1: Vec3, p2: Vec3, q2: Vec3) [2]f32 {
    const d1 = q1.sub(p1);
    const d2 = q2.sub(p2);
    const r = p1.sub(p2);
    const a = d1.dot(d1);
    const e = d2.dot(d2);
    const f = d2.dot(r);
    const eps = 1e-12;
    if (a <= eps and e <= eps) return .{ 0, 0 };
    if (a <= eps) return .{ 0, std.math.clamp(f / e, 0, 1) };
    const c = d1.dot(r);
    if (e <= eps) return .{ std.math.clamp(-c / a, 0, 1), 0 };
    const b = d1.dot(d2);
    const denom = a * e - b * b;
    var s: f32 = if (denom > eps) std.math.clamp((b * f - c * e) / denom, 0, 1) else 0;
    var t = (b * s + f) / e;
    if (t < 0) {
        t = 0;
        s = std.math.clamp(-c / a, 0, 1);
    } else if (t > 1) {
        t = 1;
        s = std.math.clamp((b - c) / a, 0, 1);
    }
    return .{ s, t };
}

/// The fraction along `a b` of its point nearest `p`.
pub fn closestOnSegment(a: Vec3, b: Vec3, p: Vec3) f32 {
    const d = b.sub(a);
    const len_sq = d.lenSq();
    if (len_sq < 1e-20) return 0;
    return std.math.clamp(p.sub(a).dot(d) / len_sq, 0, 1);
}

test "a transform and its inverse undo each other, and two compose as applied one after the other" {
    const t: Transform = .{ .position = .init(1, 2, 3), .rotation = .fromAxisAngle(.unit_y, 0.7) };
    const p = Vec3.init(-4, 5, 0.5);
    try testing.expect(t.unapply(t.apply(p)).approxEql(p));
    try testing.expect(t.inverse().apply(t.apply(p)).approxEql(p));
    const u: Transform = .{ .position = .init(0, -1, 2), .rotation = .fromAxisAngle(.unit_x, -1.2) };
    try testing.expect(t.mul(u).apply(p).approxEql(t.apply(u.apply(p))));
    try testing.expect(t.invMul(u).apply(p).approxEql(t.unapply(u.apply(p))));
}

test "a ray enters a box where it first crosses a face, or misses it" {
    const box: Aabb = .{ .min = .init(-1, -1, -1), .max = .init(1, 1, 1) };
    try testing.expectApproxEqAbs(@as(f32, 0.4), box.rayCast(.init(-5, 0, 0), .init(10, 0, 0), 1).?, 1e-6);
    try testing.expect(box.rayCast(.init(-5, 2, 0), .init(10, 0, 0), 1) == null);
    try testing.expect(box.rayCast(.init(-5, 0, 0), .init(10, 0, 0), 0.3) == null);
    try testing.expectEqual(@as(f32, 0), box.rayCast(.zero, .init(0, 3, 0), 1).?);
}

test "the closest points of two segments, crossing and parallel" {
    const st = closestSegmentSegment(.init(-1, 0, 0), .init(1, 0, 0), .init(0, -1, 1), .init(0, 1, 1));
    try testing.expectApproxEqAbs(@as(f32, 0.5), st[0], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.5), st[1], 1e-6);
    const n = basis(.init(0, 1, 0));
    try testing.expectApproxEqAbs(@as(f32, 0), n[0].dot(.init(0, 1, 0)), 1e-6);
    try testing.expect(n[0].cross(n[1]).approxEql(.init(0, 1, 0)));
}

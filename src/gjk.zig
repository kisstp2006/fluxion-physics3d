// SPDX-License-Identifier: BSD-2-Clause

//! How far apart two convex solids are, and where they are nearest: the
//! Gilbert-Johnson-Keerthi distance algorithm.
//!
//! **Each solid is a core and a radius.** A ball is a point and its radius,
//! a capsule a segment and its radius, a box or a hull its corners and no
//! radius. The algorithm works on the cores alone - the hull of a few
//! points - and the radii are taken off the distance after, which keeps it
//! exact for the round shapes and quick for all of them.
//!
//! **It walks the Minkowski difference**, the set of every point of one
//! core taken from every point of the other, without ever making it: a
//! simplex of up to four of its points, each found as the furthest point of
//! the difference along a direction - the furthest corner of one solid one
//! way and of the other the other way - is shrunk to the part nearest the
//! origin and grown towards it until it can get no nearer. The origin
//! inside it means the cores overlap.
//!
//! **A cast is a distance asked again and again** - conservative
//! advancement: at each distance the solid may move that far towards the
//! other and no collision can have been missed, so it does, until the gap
//! is a target, or it has gone the whole way.

const std = @import("std");
const testing = std.testing;

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Transform = geometry.Transform;

/// A convex core: the hull of `points` in the frame `xf` places, and how
/// far round it the solid reaches.
pub const Proxy = struct {
    points: []const Vec3,
    xf: Transform = .identity,
    radius: f32 = 0,

    /// The index of the core's point furthest along `d`, a direction
    /// outside.
    pub fn support(self: Proxy, d: Vec3) u32 {
        const local = self.xf.unturn(d);
        var best: u32 = 0;
        var most = self.points[0].dot(local);
        for (self.points[1..], 1..) |p, i| {
            const along = p.dot(local);
            if (along > most) {
                most = along;
                best = @intCast(i);
            }
        }
        return best;
    }

    pub fn point(self: Proxy, i: u32) Vec3 {
        return self.xf.apply(self.points[i]);
    }
};

const Vertex = struct {
    /// The point of each core, and the second less the first.
    a: Vec3,
    b: Vec3,
    w: Vec3,
    ia: u32,
    ib: u32,
    /// Its weight in the point of the simplex nearest the origin.
    weight: f32 = 0,
};

pub const Result = struct {
    /// The nearest points of the two cores, and how far apart they are.
    point_a: Vec3,
    point_b: Vec3,
    distance: f32,
    /// Whether the cores overlap: then the points mean nothing.
    overlap: bool,
    iterations: u32,

    /// The gap between the two solids, radii and all: negative is not a
    /// depth, only that they touch.
    pub fn separation(self: Result, a: Proxy, b: Proxy) f32 {
        return self.distance - a.radius - b.radius;
    }

    /// From `a` to `b`, where the cores are apart.
    pub fn normal(self: Result) Vec3 {
        return self.point_b.sub(self.point_a).tryNorm() orelse Vec3.init(0, 1, 0);
    }
};

const max_iterations = 48;

/// How near the cores of `a` and `b` come.
pub fn distance(a: Proxy, b: Proxy) Result {
    var simplex: [4]Vertex = undefined;
    var count: usize = 1;
    // A first vertex from the first points; the rest along the gap.
    simplex[0] = vertexOf(a, b, 0, 0);

    var iterations: u32 = 0;
    var last_distance_sq: f32 = std.math.inf(f32);
    while (iterations < max_iterations) : (iterations += 1) {
        const closest = solve(&simplex, &count);
        if (count == 4) return finish(simplex[0..count], true, iterations);
        const distance_sq = closest.lenSq();
        if (distance_sq < 1e-14) return finish(simplex[0..count], true, iterations);
        // No nearer than last time: it has converged.
        if (distance_sq >= last_distance_sq) break;
        last_distance_sq = distance_sq;

        const d = closest.neg();
        const next = vertexOf(a, b, a.support(d.neg()), b.support(d));
        // A point already in it, or one no further along: the nearest has
        // been found.
        const repeated = for (simplex[0..count]) |v| {
            if (v.ia == next.ia and v.ib == next.ib) break true;
        } else false;
        if (repeated) break;
        if (next.w.dot(d) - closest.dot(d) <= 1e-7 * @max(1, distance_sq)) {
            // Progress along the direction no greater than rounding.
            if (closest.dot(d) - next.w.dot(d) >= -1e-6 * @sqrt(distance_sq)) break;
        }
        simplex[count] = next;
        count += 1;
    }
    _ = solve(&simplex, &count);
    return finish(simplex[0..count], false, iterations);
}

fn vertexOf(a: Proxy, b: Proxy, ia: u32, ib: u32) Vertex {
    const pa = a.point(ia);
    const pb = b.point(ib);
    return .{ .a = pa, .b = pb, .w = pb.sub(pa), .ia = ia, .ib = ib };
}

fn finish(simplex: []const Vertex, overlap: bool, iterations: u32) Result {
    var pa: Vec3 = .zero;
    var pb: Vec3 = .zero;
    for (simplex) |v| {
        pa = pa.mulAdd(v.a, v.weight);
        pb = pb.mulAdd(v.b, v.weight);
    }
    return .{ .point_a = pa, .point_b = pb, .distance = if (overlap) 0 else pa.dist(pb), .overlap = overlap, .iterations = iterations };
}

/// Shrink the simplex to the part of it nearest the origin, weight each
/// point it keeps, and hand back that nearest point.
fn solve(simplex: *[4]Vertex, count: *usize) Vec3 {
    switch (count.*) {
        1 => {
            simplex[0].weight = 1;
            return simplex[0].w;
        },
        2 => return solveSegment(simplex, count, 0, 1),
        3 => return solveTriangle(simplex, count, 0, 1, 2),
        else => return solveTetrahedron(simplex, count),
    }
}

/// Keep just the points `keep` of the simplex, in that order.
fn keepOnly(simplex: *[4]Vertex, count: *usize, keep: []const usize) void {
    var kept: [4]Vertex = undefined;
    for (keep, 0..) |k, i| kept[i] = simplex[k];
    for (keep, 0..) |_, i| simplex[i] = kept[i];
    count.* = keep.len;
}

fn solveSegment(simplex: *[4]Vertex, count: *usize, i: usize, j: usize) Vec3 {
    const a = simplex[i].w;
    const b = simplex[j].w;
    const ab = b.sub(a);
    const t = -a.dot(ab);
    if (t <= 0) {
        keepOnly(simplex, count, &.{i});
        simplex[0].weight = 1;
        return simplex[0].w;
    }
    const len_sq = ab.lenSq();
    if (t >= len_sq) {
        keepOnly(simplex, count, &.{j});
        simplex[0].weight = 1;
        return simplex[0].w;
    }
    const s = t / len_sq;
    keepOnly(simplex, count, &.{ i, j });
    simplex[0].weight = 1 - s;
    simplex[1].weight = s;
    return a.mulAdd(ab, s);
}

/// Ericson's nearest point of a triangle, by the region of it the origin
/// is nearest.
fn solveTriangle(simplex: *[4]Vertex, count: *usize, i: usize, j: usize, k: usize) Vec3 {
    const a = simplex[i].w;
    const b = simplex[j].w;
    const c = simplex[k].w;
    const ab = b.sub(a);
    const ac = c.sub(a);
    const ap = a.neg();
    const d1 = ab.dot(ap);
    const d2 = ac.dot(ap);
    if (d1 <= 0 and d2 <= 0) return one(simplex, count, i);
    const bp = b.neg();
    const d3 = ab.dot(bp);
    const d4 = ac.dot(bp);
    if (d3 >= 0 and d4 <= d3) return one(simplex, count, j);
    const vc = d1 * d4 - d3 * d2;
    if (vc <= 0 and d1 >= 0 and d3 <= 0) {
        const v = d1 / (d1 - d3);
        return two(simplex, count, i, j, v);
    }
    const cp = c.neg();
    const d5 = ab.dot(cp);
    const d6 = ac.dot(cp);
    if (d6 >= 0 and d5 <= d6) return one(simplex, count, k);
    const vb = d5 * d2 - d1 * d6;
    if (vb <= 0 and d2 >= 0 and d6 <= 0) {
        const w = d2 / (d2 - d6);
        return two(simplex, count, i, k, w);
    }
    const va = d3 * d6 - d5 * d4;
    if (va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0) {
        const w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return two(simplex, count, j, k, w);
    }
    const denom = 1 / (va + vb + vc);
    const v = vb * denom;
    const w = vc * denom;
    keepOnly(simplex, count, &.{ i, j, k });
    simplex[0].weight = 1 - v - w;
    simplex[1].weight = v;
    simplex[2].weight = w;
    return a.add(ab.scale(v)).add(ac.scale(w));
}

fn one(simplex: *[4]Vertex, count: *usize, i: usize) Vec3 {
    keepOnly(simplex, count, &.{i});
    simplex[0].weight = 1;
    return simplex[0].w;
}

fn two(simplex: *[4]Vertex, count: *usize, i: usize, j: usize, t: f32) Vec3 {
    keepOnly(simplex, count, &.{ i, j });
    simplex[0].weight = 1 - t;
    simplex[1].weight = t;
    return simplex[0].w.lerp(simplex[1].w, t);
}

/// The origin inside the tetrahedron, or the nearest point of whichever of
/// its faces it is outside of that is nearest.
fn solveTetrahedron(simplex: *[4]Vertex, count: *usize) Vec3 {
    const faces = [_][4]usize{ .{ 0, 1, 2, 3 }, .{ 0, 2, 3, 1 }, .{ 0, 3, 1, 2 }, .{ 1, 3, 2, 0 } };
    var best: ?Vec3 = null;
    var best_simplex: [4]Vertex = undefined;
    var best_count: usize = 0;
    for (faces) |f| {
        const a = simplex[f[0]].w;
        const b = simplex[f[1]].w;
        const c = simplex[f[2]].w;
        const opposite = simplex[f[3]].w;
        const n = b.sub(a).cross(c.sub(a));
        // The origin and the fourth point on the same side: not this face.
        const origin_side = n.dot(a.neg());
        const other_side = n.dot(opposite.sub(a));
        if (origin_side * other_side > 0) continue;
        var trial = simplex.*;
        var trial_count: usize = 4;
        const near = solveTriangle(&trial, &trial_count, f[0], f[1], f[2]);
        if (best == null or near.lenSq() < best.?.lenSq()) {
            best = near;
            best_simplex = trial;
            best_count = trial_count;
        }
    }
    const found = best orelse {
        // Inside every face: the cores overlap.
        count.* = 4;
        return .zero;
    };
    simplex.* = best_simplex;
    count.* = best_count;
    return found;
}

pub const CastResult = struct {
    /// How far along the translation it first comes within the target of
    /// touching.
    fraction: f32,
    /// Where, and out of `a` towards `b`, there.
    point: Vec3,
    normal: Vec3,
    /// Overlapping at the start: no direction to tell.
    initially_overlapping: bool,
};

/// Where `b`, moved along `translation`, first comes within `target` of
/// touching `a` - no further than `max_fraction` of the way - or null.
pub fn cast(a: Proxy, b: Proxy, translation: Vec3, target: f32, max_fraction: f32) ?CastResult {
    var moved = b;
    var t: f32 = 0;
    const radii = a.radius + b.radius;
    var iterations: u32 = 0;
    while (iterations < 40) : (iterations += 1) {
        moved.xf.position = b.xf.position.mulAdd(translation, t);
        const found = distance(a, moved);
        if (found.overlap) {
            if (t == 0) return .{ .fraction = 0, .point = found.point_a, .normal = translation.neg().tryNorm() orelse Vec3.init(0, 1, 0), .initially_overlapping = true };
            return null;
        }
        const gap = found.distance - radii;
        const n = found.normal();
        if (gap <= target) {
            if (t == 0 and gap < 0) return .{ .fraction = 0, .point = found.point_a.mulAdd(n, a.radius), .normal = n, .initially_overlapping = true };
            return .{ .fraction = t, .point = found.point_a.mulAdd(n, a.radius), .normal = n, .initially_overlapping = false };
        }
        // How fast the gap closes along the line between the nearest points.
        const closing = -translation.dot(n);
        if (closing <= 0) return null;
        t += (gap - target * 0.5) / closing;
        if (t > max_fraction) return null;
    }
    return null;
}

const cube = [_]Vec3{
    .init(-1, -1, -1), .init(1, -1, -1), .init(-1, 1, -1), .init(1, 1, -1),
    .init(-1, -1, 1),  .init(1, -1, 1),  .init(-1, 1, 1),  .init(1, 1, 1),
};

test "the distance between two boxes apart, and two that overlap" {
    const a: Proxy = .{ .points = &cube };
    var b: Proxy = .{ .points = &cube, .xf = .at(.init(5, 0.5, 0)) };
    const apart = distance(a, b);
    try testing.expect(!apart.overlap);
    try testing.expectApproxEqAbs(@as(f32, 3), apart.distance, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, 1), apart.point_a.x, 1e-5);
    try testing.expect(apart.normal().approxEql(.init(1, 0, 0)));
    b.xf.position = .init(1.5, 0.3, 0.2);
    try testing.expect(distance(a, b).overlap);
    // Turned an eighth about z: its edge towards the first's face.
    b.xf = .{ .position = .init(4, 0, 0), .rotation = .fromAxisAngle(.init(0, 0, 1), std.math.pi / 4.0) };
    try testing.expectApproxEqAbs(@as(f32, 3 - std.math.sqrt2), distance(a, b).distance, 1e-4);
}

test "a ball and a capsule are a point and a segment with their radii" {
    const centre = [_]Vec3{.zero};
    const segment = [_]Vec3{ .init(0, -1, 0), .init(0, 1, 0) };
    const ball: Proxy = .{ .points = &centre, .xf = .at(.init(3, 0.5, 0)), .radius = 0.5 };
    const capsule: Proxy = .{ .points = &segment, .radius = 0.5 };
    const found = distance(capsule, ball);
    try testing.expectApproxEqAbs(@as(f32, 3), found.distance, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, 2), found.separation(capsule, ball), 1e-5);
    try testing.expect(found.point_a.approxEql(.init(0, 0.5, 0)));
}

test "a box cast at another stops a target short of it, and one cast away misses" {
    const a: Proxy = .{ .points = &cube };
    const b: Proxy = .{ .points = &cube, .xf = .at(.init(10, 0, 0)) };
    const hit = cast(a, b, .init(-20, 0, 0), 0.01, 1).?;
    // From ten to two apart, less the target: eight of twenty.
    try testing.expectApproxEqAbs(@as(f32, (8 - 0.01) / 20.0), hit.fraction, 1e-3);
    try testing.expect(hit.normal.approxEql(.init(1, 0, 0)));
    try testing.expect(cast(a, b, .init(20, 0, 0), 0.01, 1) == null);
    try testing.expect(cast(a, b, .init(-5, 0, 0), 0.01, 1) == null);
}

// SPDX-License-Identifier: BSD-2-Clause

//! Where two solids touch, and how deep: a normal and up to four points,
//! each with how far apart the two surfaces are there - negative where they
//! have sunk in. Points are kept while the surfaces are within a margin of
//! each other, not only when they have met, so the solver can stop a fast
//! body at the surface rather than finding it already through.
//!
//! **Round solids are cores with radii.** A ball against a ball, a capsule
//! against a capsule, is the nearest points of two points or segments; two
//! capsules lying side by side touch along a line, at both ends of it.
//!
//! **A round solid against one with corners** is the nearest points of its
//! core and the solid, by `gjk` - and a capsule lying along a face touches
//! it at the two ends of the part of it over the face. A core that has gone
//! inside is pushed out the shortest way, found by trying each face and,
//! for a capsule, each edge crossed with it.
//!
//! **Two solids with corners** are tried along every way they could be
//! apart - each face of either, and each edge of one crossed with each
//! edge of the other - and the way they are least sunk along is the one
//! they are pushed apart along: Gottschalk's separating axis test, made
//! quick for hulls by Gregorius's rule that only edges whose faces' normals
//! cross on the sphere of directions can be nearest. Along a face, the face
//! of the other solid that faces it most is clipped to it - Sutherland and
//! Hodgman - and its corners that stay are the points; four, the deepest
//! and the three that span the most, are kept. Along two edges, it is the
//! nearest points of the two.
//!
//! Each point has an `id` from the faces and corners it came from, so the
//! same point the next step is known, and starts from the push it ended
//! with.

const std = @import("std");
const testing = std.testing;

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Transform = geometry.Transform;
const hull_mod = @import("hull.zig");
const Poly = hull_mod.Poly;
const Face = hull_mod.Face;
const shape_mod = @import("shape.zig");
const Geometry = shape_mod.Geometry;
const gjk = @import("gjk.zig");

pub const max_points = 4;

pub const ManifoldPoint = struct {
    /// In the world, half way between the two surfaces.
    point: Vec3,
    /// Apart, or sunk in where negative.
    separation: f32,
    id: u32,
};

pub const Manifold = struct {
    /// From the first solid to the second.
    normal: Vec3 = .zero,
    points: [max_points]ManifoldPoint = undefined,
    count: u32 = 0,

    pub fn pointSlice(self: *const Manifold) []const ManifoldPoint {
        return self.points[0..self.count];
    }

    fn add(self: *Manifold, p: ManifoldPoint) void {
        if (self.count == max_points) return;
        self.points[self.count] = p;
        self.count += 1;
    }

    /// The least apart its points are.
    pub fn deepest(self: *const Manifold) f32 {
        var out: f32 = std.math.inf(f32);
        for (self.pointSlice()) |p| out = @min(out, p.separation);
        return out;
    }
};

/// A solid placed in the world, as collide sees it.
pub const Solid = union(enum) {
    round: Round,
    poly: Placed,
};

/// A point or a segment, and how far round it.
pub const Round = struct {
    core: [2]Vec3,
    count: u32,
    radius: f32,

    fn proxy(self: *const Round) gjk.Proxy {
        return .{ .points = self.core[0..self.count] };
    }
};

/// A solid with corners, and where it is.
pub const Placed = struct {
    poly: Poly,
    xf: Transform,

    fn vertex(self: Placed, i: usize) Vec3 {
        return self.xf.apply(self.poly.vertices[i]);
    }

    fn normalOf(self: Placed, f: usize) Vec3 {
        return self.xf.turn(self.poly.faces[f].normal);
    }

    /// Face `f`'s plane in the world: `normal . x = distance`.
    fn plane(self: Placed, f: usize) struct { normal: Vec3, distance: f32 } {
        const n = self.normalOf(f);
        return .{ .normal = n, .distance = self.poly.faces[f].distance + n.dot(self.xf.position) };
    }

    /// Its corner furthest along `d`, in the world.
    fn support(self: Placed, d: Vec3) Vec3 {
        return self.vertex(self.poly.support(self.xf.unturn(d)));
    }

    fn center(self: Placed) Vec3 {
        var sum: Vec3 = .zero;
        for (self.poly.vertices) |v| sum = sum.add(v);
        return self.xf.apply(sum.scale(1 / @as(f32, @floatFromInt(self.poly.vertices.len))));
    }

    fn proxy(self: Placed) gjk.Proxy {
        return .{ .points = self.poly.vertices, .xf = self.xf };
    }

    /// The face whose normal is most along `n`.
    fn faceAlong(self: Placed, n: Vec3) usize {
        const local = self.xf.unturn(n);
        var best: usize = 0;
        var most: f32 = -std.math.inf(f32);
        for (self.poly.faces, 0..) |f, i| {
            const along = f.normal.dot(local);
            if (along > most) {
                most = along;
                best = i;
            }
        }
        return best;
    }

    /// How big it is: the length an edge's line is stretched to.
    fn reach(self: Placed) f32 {
        var most: f32 = 0;
        for (self.poly.vertices) |v| most = @max(most, v.lenSq());
        return 2 * @sqrt(most) + 1;
    }
};

/// Room for the solids `solidOf` makes of a box, a cylinder or a triangle,
/// and the cores `proxyOf` makes of a ball and a capsule.
pub const Storage = struct {
    box: hull_mod.BoxPoly = undefined,
    cylinder: hull_mod.CylinderPoly = undefined,
    triangle: hull_mod.TrianglePoly = undefined,
    core: [3]Vec3 = undefined,
};

/// A geometry placed by `xf` as a core and a radius, for `gjk` - or null
/// for a mesh.
pub fn proxyOf(g: Geometry, xf: Transform, storage: *Storage) ?gjk.Proxy {
    switch (g) {
        .sphere => |s| {
            storage.core[0] = .zero;
            return .{ .points = storage.core[0..1], .xf = xf, .radius = s.radius };
        },
        .capsule => |c| {
            storage.core[0] = .init(0, -c.half_height, 0);
            storage.core[1] = .init(0, c.half_height, 0);
            return .{ .points = storage.core[0..2], .xf = xf, .radius = c.radius };
        },
        .box => |b| {
            storage.box = .init(b.half);
            return .{ .points = &storage.box.vertices, .xf = xf };
        },
        .cylinder => |c| {
            storage.cylinder = .init(c.half_height, c.radius);
            return .{ .points = &storage.cylinder.vertices, .xf = xf };
        },
        .hull => |h| return .{ .points = h.vertices, .xf = xf },
        .mesh => return null,
    }
}

/// Triangle `a b c`, in the world, as a core.
pub fn triangleProxy(storage: *Storage, a: Vec3, b: Vec3, c: Vec3) gjk.Proxy {
    storage.core = .{ a, b, c };
    return .{ .points = &storage.core };
}

/// A geometry placed by `xf`, as a solid - or null for a mesh, which is
/// its triangles: see `triangleSolid`.
pub fn solidOf(g: Geometry, xf: Transform, storage: *Storage) ?Solid {
    switch (g) {
        .sphere => |s| return .{ .round = .{ .core = .{ xf.position, xf.position }, .count = 1, .radius = s.radius } },
        .capsule => |c| {
            const axis = xf.turn(.init(0, c.half_height, 0));
            return .{ .round = .{ .core = .{ xf.position.sub(axis), xf.position.add(axis) }, .count = 2, .radius = c.radius } };
        },
        .box => |b| {
            storage.box = .init(b.half);
            return .{ .poly = .{ .poly = storage.box.poly(), .xf = xf } };
        },
        .cylinder => |c| {
            storage.cylinder = .init(c.half_height, c.radius);
            return .{ .poly = .{ .poly = storage.cylinder.poly(), .xf = xf } };
        },
        .hull => |h| return .{ .poly = .{ .poly = h.poly(), .xf = xf } },
        .mesh => return null,
    }
}

/// Triangle `a b c`, in the world, as a solid.
pub fn triangleSolid(storage: *Storage, a: Vec3, b: Vec3, c: Vec3) Solid {
    storage.triangle = .init(a, b, c);
    return .{ .poly = .{ .poly = storage.triangle.poly(), .xf = .identity } };
}

/// Where `a` and `b` touch, or come within `margin` of it, the normal from
/// `a` to `b`. No points: further apart than that.
pub fn collide(a: Solid, b: Solid, margin: f32) Manifold {
    return switch (a) {
        .round => |ra| switch (b) {
            .round => |rb| roundRound(ra, rb, margin),
            .poly => |pb| flipped(roundPoly(ra, pb, margin)),
        },
        .poly => |pa| switch (b) {
            .round => |rb| roundPoly(rb, pa, margin),
            .poly => |pb| polyPoly(pa, pb, margin),
        },
    };
}

fn flipped(m: Manifold) Manifold {
    var out = m;
    out.normal = m.normal.neg();
    return out;
}

/// Any direction at right angles to `v`, or up.
fn across(v: Vec3) Vec3 {
    return (v.anyPerp().tryNorm()) orelse Vec3.init(0, 1, 0);
}

// -------------------------------------------------------------------------
// Round against round
// -------------------------------------------------------------------------

fn roundRound(a: Round, b: Round, margin: f32) Manifold {
    var m: Manifold = .{};
    const radii = a.radius + b.radius;

    // Two segments nearly parallel: the two ends of their overlap.
    if (a.count == 2 and b.count == 2) parallel: {
        const da = a.core[1].sub(a.core[0]);
        const db = b.core[1].sub(b.core[0]);
        const la = da.len();
        const lb = db.len();
        if (la < 1e-6 or lb < 1e-6) break :parallel;
        const ua = da.scale(1 / la);
        if (@abs(ua.dot(db.scale(1 / lb))) < 0.995) break :parallel;
        const t0 = std.math.clamp(b.core[0].sub(a.core[0]).dot(ua), 0, la);
        const t1 = std.math.clamp(b.core[1].sub(a.core[0]).dot(ua), 0, la);
        const lo = @min(t0, t1);
        const hi = @max(t0, t1);
        if (hi - lo < 1e-3 * la) break :parallel;
        // One normal for both: from A's line to B's, at the middle.
        const mid_a = a.core[0].mulAdd(ua, (lo + hi) / 2);
        const mid_b = b.core[0].lerp(b.core[1], geometry.closestOnSegment(b.core[0], b.core[1], mid_a));
        const gap = mid_b.sub(mid_a);
        const n = if (gap.lenSq() > 1e-12) gap.norm() else across(ua);
        m.normal = n;
        for ([_]f32{ lo, hi }, 0..) |t, i| {
            const pa = a.core[0].mulAdd(ua, t);
            const pb = b.core[0].lerp(b.core[1], geometry.closestOnSegment(b.core[0], b.core[1], pa));
            const sep = pb.sub(pa).dot(n) - radii;
            if (sep > margin) continue;
            m.add(.{ .point = pa.mulAdd(n, a.radius + sep / 2), .separation = sep, .id = @intCast(i) });
        }
        return m;
    }

    const pa, const pb = nearest(a, b);
    const gap = pb.sub(pa);
    const dist = gap.len();
    const sep = dist - radii;
    if (sep > margin) return m;
    const n = if (dist > 1e-6) gap.scale(1 / dist) else if (a.count == 2) across(a.core[1].sub(a.core[0])) else Vec3.init(0, 1, 0);
    m.normal = n;
    m.add(.{ .point = pa.mulAdd(n, a.radius + sep / 2), .separation = sep, .id = 0 });
    return m;
}

/// The nearest points of two cores.
fn nearest(a: Round, b: Round) [2]Vec3 {
    if (a.count == 1 and b.count == 1) return .{ a.core[0], b.core[0] };
    if (a.count == 1) return .{ a.core[0], b.core[0].lerp(b.core[1], geometry.closestOnSegment(b.core[0], b.core[1], a.core[0])) };
    if (b.count == 1) return .{ a.core[0].lerp(a.core[1], geometry.closestOnSegment(a.core[0], a.core[1], b.core[0])), b.core[0] };
    const st = geometry.closestSegmentSegment(a.core[0], a.core[1], b.core[0], b.core[1]);
    return .{ a.core[0].lerp(a.core[1], st[0]), b.core[0].lerp(b.core[1], st[1]) };
}

// -------------------------------------------------------------------------
// Round against a solid with corners: the normal out of the solid
// -------------------------------------------------------------------------

fn roundPoly(r: Round, p: Placed, margin: f32) Manifold {
    var m: Manifold = .{};
    const found = gjk.distance(p.proxy(), r.proxy());
    if (!found.overlap and found.distance > 1e-5) {
        const sep = found.distance - r.radius;
        if (sep > margin) return m;
        const n = found.point_b.sub(found.point_a).scale(1 / found.distance);
        // A segment lying along a face touches it at both ends.
        if (r.count == 2) {
            const face = p.faceAlong(n);
            const plane = p.plane(face);
            const length = r.core[1].sub(r.core[0]).len();
            if (plane.normal.dot(n) > 0.98 and @abs(plane.normal.dot(r.core[1].sub(r.core[0]))) < 0.15 * length) {
                if (segmentOnFace(r, p, face, margin, &m)) return m;
            }
        }
        m.normal = n;
        m.add(.{ .point = found.point_a.add(found.point_b.mulAdd(n, -r.radius)).scale(0.5), .separation = sep, .id = 0 });
        return m;
    }
    return roundInside(r, p, margin);
}

/// The part of a segment over face `face`, its two ends each a point, the
/// normal the face's. False when none of it is over the face.
fn segmentOnFace(r: Round, p: Placed, face: usize, margin: f32, m: *Manifold) bool {
    const plane = p.plane(face);
    const corners = p.poly.cornersOf(p.poly.faces[face]);
    var c0 = r.core[0];
    var c1 = r.core[1];
    for (corners, 0..) |v, k| {
        const from = p.vertex(v);
        const to = p.vertex(corners[(k + 1) % corners.len]);
        const side = to.sub(from).cross(plane.normal);
        const offset = side.dot(from);
        const d0 = side.dot(c0) - offset;
        const d1 = side.dot(c1) - offset;
        if (d0 > 0 and d1 > 0) return false;
        if (d0 > 0) c0 = c0.lerp(c1, d0 / (d0 - d1));
        if (d1 > 0) c1 = c1.lerp(c0, d1 / (d1 - d0));
    }
    m.normal = plane.normal;
    for ([_]Vec3{ c0, c1 }, 0..) |q, i| {
        const sep = plane.normal.dot(q) - plane.distance - r.radius;
        if (sep > margin) continue;
        m.add(.{ .point = q.mulAdd(plane.normal, -(r.radius + sep / 2)), .separation = sep, .id = @intCast(i) });
    }
    return m.count > 0;
}

/// A core sunk inside a solid: out through whichever face - or for a
/// segment, across whichever edge - is nearest.
fn roundInside(r: Round, p: Placed, margin: f32) Manifold {
    var m: Manifold = .{};
    const core = r.core[0..r.count];
    var best_sep: f32 = -std.math.inf(f32);
    var best_face: usize = 0;
    for (p.poly.faces, 0..) |_, f| {
        const plane = p.plane(f);
        var low: f32 = std.math.inf(f32);
        for (core) |c| low = @min(low, plane.normal.dot(c) - plane.distance);
        if (low > best_sep) {
            best_sep = low;
            best_face = f;
        }
    }
    var edge_axis: ?Vec3 = null;
    var edge_support: Vec3 = .zero;
    var edge_dir: Vec3 = .zero;
    if (r.count == 2) {
        const u = core[1].sub(core[0]);
        const middle = core[0].add(core[1]).scale(0.5);
        const center = p.center();
        for (p.poly.edges) |e| {
            const dir = p.vertex(e[1]).sub(p.vertex(e[0]));
            var axis = dir.cross(u).tryNorm() orelse continue;
            if (axis.dot(middle.sub(center)) < 0) axis = axis.neg();
            const far = p.support(axis);
            const sep = @min(axis.dot(core[0]), axis.dot(core[1])) - axis.dot(far);
            // Faces win ties: their points are steadier.
            if (sep > best_sep + 1e-4) {
                best_sep = sep;
                edge_axis = axis;
                edge_support = far;
                edge_dir = dir;
            }
        }
    }
    if (best_sep - r.radius > margin) return m;

    if (edge_axis) |axis| {
        const reach = p.reach();
        const unit = edge_dir.norm();
        const st = geometry.closestSegmentSegment(core[0], core[1], edge_support.mulAdd(unit, -reach), edge_support.mulAdd(unit, reach));
        const on_core = core[0].lerp(core[1], st[0]);
        const sep = best_sep - r.radius;
        m.normal = axis;
        m.add(.{ .point = on_core.mulAdd(axis, -(r.radius + sep / 2)), .separation = sep, .id = 2 });
        return m;
    }
    if (r.count == 2 and segmentOnFace(r, p, best_face, margin, &m)) return m;
    const plane = p.plane(best_face);
    // The core's deepest point.
    var deepest = core[0];
    for (core) |c| if (plane.normal.dot(c) < plane.normal.dot(deepest)) {
        deepest = c;
    };
    const sep = plane.normal.dot(deepest) - plane.distance - r.radius;
    m.normal = plane.normal;
    m.add(.{ .point = deepest.mulAdd(plane.normal, -(r.radius + sep / 2)), .separation = sep, .id = 0 });
    return m;
}

// -------------------------------------------------------------------------
// Two solids with corners
// -------------------------------------------------------------------------

const FaceQuery = struct { separation: f32, face: usize };

/// The face of `a` that `b` is furthest out of.
fn faceQuery(a: Placed, b: Placed) FaceQuery {
    var best: FaceQuery = .{ .separation = -std.math.inf(f32), .face = 0 };
    for (a.poly.faces, 0..) |_, f| {
        const plane = a.plane(f);
        const sep = plane.normal.dot(b.support(plane.normal.neg())) - plane.distance;
        if (sep > best.separation) best = .{ .separation = sep, .face = f };
    }
    return best;
}

const EdgeQuery = struct {
    separation: f32,
    /// Out of `a` towards `b`.
    axis: Vec3,
    /// A point on each edge, and which way it runs.
    point_a: Vec3,
    dir_a: Vec3,
    point_b: Vec3,
    dir_b: Vec3,
    id: u32,
};

/// The pair of edges, one of each, that the two are furthest apart across.
fn edgeQuery(a: Placed, b: Placed) EdgeQuery {
    var best: EdgeQuery = .{ .separation = -std.math.inf(f32), .axis = .zero, .point_a = .zero, .dir_a = .zero, .point_b = .zero, .dir_b = .zero, .id = 0 };
    const center_a = a.center();
    const center_b = b.center();
    const pruned = a.poly.edge_faces.len == a.poly.edges.len and b.poly.edge_faces.len == b.poly.edges.len;
    for (a.poly.edges, 0..) |ea, i| {
        const pa = a.vertex(ea[0]);
        const da = a.vertex(ea[1]).sub(pa);
        for (b.poly.edges, 0..) |eb, j| {
            const pb = b.vertex(eb[0]);
            const db = b.vertex(eb[1]).sub(pb);
            var axis = da.cross(db).tryNorm() orelse continue;
            var sep: f32 = undefined;
            var on_a = pa;
            var on_b = pb;
            if (pruned) {
                // Only edges whose faces' normals cross on the sphere of
                // directions can be nearest: then the edges themselves are.
                const fa = a.poly.edge_faces[i];
                const fb = b.poly.edge_faces[j];
                if (!isMinkowskiFace(a.normalOf(fa[0]), a.normalOf(fa[1]), b.normalOf(fb[0]).neg(), b.normalOf(fb[1]).neg())) continue;
                if (axis.dot(pa.sub(center_a)) < 0) axis = axis.neg();
                sep = axis.dot(pb.sub(pa));
            } else {
                // Edges standing for every edge their way: how far apart the
                // two are along the axis, whichever edges those are.
                if (axis.dot(center_b.sub(center_a)) < 0) axis = axis.neg();
                on_a = a.support(axis);
                on_b = b.support(axis.neg());
                sep = axis.dot(on_b) - axis.dot(on_a);
            }
            if (sep > best.separation) best = .{
                .separation = sep,
                .axis = axis,
                .point_a = on_a,
                .dir_a = da,
                .point_b = on_b,
                .dir_b = db,
                .id = @intCast((i << 8) | j),
            };
        }
    }
    return best;
}

/// Whether the arc between `a` and `b` and the arc between `c` and `d`
/// cross on the sphere of directions.
fn isMinkowskiFace(a: Vec3, b: Vec3, c: Vec3, d: Vec3) bool {
    const bxa = b.cross(a);
    const dxc = d.cross(c);
    const cba = c.dot(bxa);
    const dba = d.dot(bxa);
    const adc = a.dot(dxc);
    const bdc = b.dot(dxc);
    return cba * dba < 0 and adc * bdc < 0 and cba * bdc > 0;
}

fn polyPoly(a: Placed, b: Placed, margin: f32) Manifold {
    const fa = faceQuery(a, b);
    if (fa.separation > margin) return .{};
    const fb = faceQuery(b, a);
    if (fb.separation > margin) return .{};
    const edge = edgeQuery(a, b);
    if (edge.separation > margin) return .{};

    // A face is the steadier contact: an edge only when it is clearly the
    // way they are least sunk.
    const tolerance = 0.002;
    const face_best = @max(fa.separation, fb.separation);
    if (edge.separation > face_best + tolerance) {
        const reach = @max(a.reach(), b.reach());
        const ua = edge.dir_a.norm();
        const ub = edge.dir_b.norm();
        const st = geometry.closestSegmentSegment(edge.point_a.mulAdd(ua, -reach), edge.point_a.mulAdd(ua, reach), edge.point_b.mulAdd(ub, -reach), edge.point_b.mulAdd(ub, reach));
        const on_a = edge.point_a.mulAdd(ua, -reach).lerp(edge.point_a.mulAdd(ua, reach), st[0]);
        const on_b = edge.point_b.mulAdd(ub, -reach).lerp(edge.point_b.mulAdd(ub, reach), st[1]);
        var m: Manifold = .{ .normal = edge.axis };
        m.add(.{ .point = on_a.add(on_b).scale(0.5), .separation = edge.separation, .id = 0x4000_0000 | edge.id });
        return m;
    }
    if (fb.separation > fa.separation + tolerance) return faceContact(b, fb.face, a, true, margin);
    return faceContact(a, fa.face, b, false, margin);
}

const ClipPoint = struct { p: Vec3, id: u32 };

/// The most corners a face is clipped with; a hull face with more is cut
/// short, which drops only points `reduce` would.
const max_clip = 64;

/// Face `ref_face` of `ref` against the face of `inc` that faces it most:
/// that face's corners clipped to it. `flip` when `ref` is the second
/// solid, so the normal goes the other way.
fn faceContact(ref: Placed, ref_face: usize, inc: Placed, flip: bool, margin: f32) Manifold {
    const plane = ref.plane(ref_face);
    const n = plane.normal;

    var inc_face: usize = 0;
    var lowest: f32 = std.math.inf(f32);
    for (inc.poly.faces, 0..) |_, f| {
        const along = inc.normalOf(f).dot(n);
        if (along < lowest) {
            lowest = along;
            inc_face = f;
        }
    }

    var buffers: [2][max_clip]ClipPoint = undefined;
    var count: usize = 0;
    for (inc.poly.cornersOf(inc.poly.faces[inc_face])) |v| {
        if (count == max_clip) break;
        buffers[0][count] = .{ .p = inc.vertex(v), .id = @intCast(count) };
        count += 1;
    }
    var current: usize = 0;
    const corners = ref.poly.cornersOf(ref.poly.faces[ref_face]);
    for (corners, 0..) |v, k| {
        const from = ref.vertex(v);
        const to = ref.vertex(corners[(k + 1) % corners.len]);
        const side = to.sub(from).cross(n);
        const offset = side.dot(from);
        count = clip(buffers[current][0..count], &buffers[1 - current], side, offset, @intCast(k));
        current = 1 - current;
        if (count == 0) break;
    }

    var candidates: [max_clip]ManifoldPoint = undefined;
    var kept: usize = 0;
    const base_id: u32 = (@as(u32, @intFromBool(flip)) << 31) | (@as(u32, @intCast(ref_face & 0xFF)) << 23) | (@as(u32, @intCast(inc_face & 0xFF)) << 15);
    for (buffers[current][0..count]) |c| {
        const sep = n.dot(c.p) - plane.distance;
        if (sep > margin) continue;
        candidates[kept] = .{ .point = c.p.mulAdd(n, -sep / 2), .separation = sep, .id = base_id | (c.id & 0x7FFF) };
        kept += 1;
    }
    var m: Manifold = .{ .normal = if (flip) n.neg() else n };
    if (kept == 0) {
        // The faces do not overlap seen along the normal: the corner of the
        // other that reaches furthest in.
        const deep = inc.support(n.neg());
        const sep = n.dot(deep) - plane.distance;
        if (sep <= margin) m.add(.{ .point = deep.mulAdd(n, -sep / 2), .separation = sep, .id = base_id | 0x7FFF });
        return m;
    }
    reduce(candidates[0..kept], n, &m);
    return m;
}

/// The part of polygon `in` on the inner side of the plane `side . x =
/// offset`, into `out`. A point made where an edge crosses the plane is
/// named by the plane and the edge.
fn clip(in: []const ClipPoint, out: *[max_clip]ClipPoint, side: Vec3, offset: f32, plane_id: u32) usize {
    var count: usize = 0;
    if (in.len == 0) return 0;
    var prev = in[in.len - 1];
    var prev_d = side.dot(prev.p) - offset;
    for (in) |cur| {
        const cur_d = side.dot(cur.p) - offset;
        if ((prev_d <= 0) != (cur_d <= 0)) {
            if (count < max_clip) {
                const t = prev_d / (prev_d - cur_d);
                out[count] = .{ .p = prev.p.lerp(cur.p, t), .id = 0x4000 | (plane_id << 7) | (prev.id & 0x7F) };
                count += 1;
            }
        }
        if (cur_d <= 0 and count < max_clip) {
            out[count] = cur;
            count += 1;
        }
        prev = cur;
        prev_d = cur_d;
    }
    return count;
}

/// Four of `points`: the deepest, the furthest from it, the one making the
/// biggest triangle with those, and the one furthest outside that.
fn reduce(points: []const ManifoldPoint, n: Vec3, m: *Manifold) void {
    if (points.len <= max_points) {
        for (points) |p| m.add(p);
        return;
    }
    var chosen: [max_points]usize = undefined;
    chosen[0] = 0;
    for (points, 0..) |p, i| if (p.separation < points[chosen[0]].separation) {
        chosen[0] = i;
    };
    const p0 = points[chosen[0]].point;
    chosen[1] = chosen[0];
    var far: f32 = -1;
    for (points, 0..) |p, i| {
        const d = p.point.distSq(p0);
        if (d > far) {
            far = d;
            chosen[1] = i;
        }
    }
    const p1 = points[chosen[1]].point;
    chosen[2] = chosen[0];
    var area: f32 = -1;
    for (points, 0..) |p, i| {
        const a = @abs(n.dot(p1.sub(p0).cross(p.point.sub(p0))));
        if (a > area) {
            area = a;
            chosen[2] = i;
        }
    }
    const p2 = points[chosen[2]].point;
    const sign: f32 = if (n.dot(p1.sub(p0).cross(p2.sub(p0))) >= 0) 1 else -1;
    var outside: f32 = 0;
    var fourth: ?usize = null;
    const tri = [_]Vec3{ p0, p1, p2 };
    for (points, 0..) |p, i| {
        if (i == chosen[0] or i == chosen[1] or i == chosen[2]) continue;
        for (0..3) |k| {
            const out = -sign * n.dot(tri[(k + 1) % 3].sub(tri[k]).cross(p.point.sub(tri[k])));
            if (out > outside) {
                outside = out;
                fourth = i;
            }
        }
    }
    m.add(points[chosen[0]]);
    if (chosen[1] != chosen[0]) m.add(points[chosen[1]]);
    if (chosen[2] != chosen[0] and chosen[2] != chosen[1]) m.add(points[chosen[2]]);
    if (fourth) |i| m.add(points[i]);
}

// -------------------------------------------------------------------------
// Tests
// -------------------------------------------------------------------------

fn solid(g: Geometry, xf: Transform, storage: *Storage) Solid {
    return solidOf(g, xf, storage).?;
}

test "two balls touch along the line between them, and not past the margin" {
    var sa: Storage = .{};
    var sb: Storage = .{};
    const a = solid(.{ .sphere = .{ .radius = 1 } }, .at(.zero), &sa);
    const b = solid(.{ .sphere = .{ .radius = 1 } }, .at(.init(1.5, 0, 0)), &sb);
    const m = collide(a, b, 0.01);
    try testing.expectEqual(@as(u32, 1), m.count);
    try testing.expect(m.normal.approxEql(.init(1, 0, 0)));
    try testing.expectApproxEqAbs(@as(f32, -0.5), m.points[0].separation, 1e-6);
    try testing.expect(m.points[0].point.approxEql(.init(0.75, 0, 0)));
    const far = solid(.{ .sphere = .{ .radius = 1 } }, .at(.init(2.1, 0, 0)), &sb);
    try testing.expectEqual(@as(u32, 0), collide(a, far, 0.05).count);
    try testing.expectEqual(@as(u32, 1), collide(a, far, 0.2).count);
}

test "a box resting on a box touches it at four corners, sunk as deep as it is" {
    var sa: Storage = .{};
    var sb: Storage = .{};
    const ground = solid(.{ .box = .{ .half = .init(5, 0.5, 5) } }, .at(.zero), &sa);
    const crate = solid(.{ .box = .{ .half = .init(0.5, 0.5, 0.5) } }, .at(.init(0.3, 0.99, -0.2)), &sb);
    const m = collide(ground, crate, 0.02);
    try testing.expectEqual(@as(u32, 4), m.count);
    try testing.expect(m.normal.approxEql(.init(0, 1, 0)));
    for (m.pointSlice()) |p| try testing.expectApproxEqAbs(@as(f32, -0.01), p.separation, 1e-5);
    // The other way round, the normal goes down.
    const back = collide(crate, ground, 0.02);
    try testing.expectEqual(@as(u32, 4), back.count);
    try testing.expect(back.normal.approxEql(.init(0, -1, 0)));
}

test "a box turned on its edge touches along that edge, and two boxes edge to edge at one point" {
    var sa: Storage = .{};
    var sb: Storage = .{};
    const ground = solid(.{ .box = .{ .half = .init(5, 0.5, 5) } }, .at(.zero), &sa);
    const on_edge = solid(.{ .box = .{ .half = .init(0.5, 0.5, 0.5) } }, .{ .position = .init(0, 0.5 + std.math.sqrt1_2 - 0.01, 0), .rotation = .fromAxisAngle(.init(0, 0, 1), std.math.pi / 4.0) }, &sb);
    const m = collide(ground, on_edge, 0.02);
    try testing.expectEqual(@as(u32, 2), m.count);
    try testing.expect(m.normal.approxEql(.init(0, 1, 0)));
    // Crossed edges: one turned about z, the other about x, edge on edge.
    var sc: Storage = .{};
    const top = solid(.{ .box = .{ .half = .init(0.5, 0.5, 0.5) } }, .{ .position = .init(0, 2 * std.math.sqrt1_2 - 0.01, 0), .rotation = .fromAxisAngle(.init(1, 0, 0), std.math.pi / 4.0) }, &sc);
    var sd: Storage = .{};
    const bottom = solid(.{ .box = .{ .half = .init(0.5, 0.5, 0.5) } }, .{ .rotation = .fromAxisAngle(.init(0, 0, 1), std.math.pi / 4.0) }, &sd);
    const crossed = collide(bottom, top, 0.02);
    try testing.expectEqual(@as(u32, 1), crossed.count);
    try testing.expect(crossed.normal.approxEql(.init(0, 1, 0)));
    try testing.expectApproxEqAbs(@as(f32, -0.01), crossed.points[0].separation, 1e-4);
}

test "a capsule lying on a box touches at both ends, and a ball sunk into one is pushed out the nearest face" {
    var sa: Storage = .{};
    var sb: Storage = .{};
    const ground = solid(.{ .box = .{ .half = .init(5, 0.5, 5) } }, .at(.zero), &sa);
    const lying = solid(.{ .capsule = .{ .half_height = 1, .radius = 0.25 } }, .{ .position = .init(0, 0.74, 0), .rotation = .fromAxisAngle(.init(0, 0, 1), std.math.pi / 2.0) }, &sb);
    const m = collide(lying, ground, 0.02);
    try testing.expectEqual(@as(u32, 2), m.count);
    try testing.expect(m.normal.approxEql(.init(0, -1, 0)));
    for (m.pointSlice()) |p| try testing.expectApproxEqAbs(@as(f32, -0.01), p.separation, 1e-5);

    const sunk = solid(.{ .sphere = .{ .radius = 0.25 } }, .at(.init(4.9, 0.2, 0)), &sb);
    const out = collide(ground, sunk, 0.02);
    try testing.expectEqual(@as(u32, 1), out.count);
    // Nearer the side than the top: out of the side.
    try testing.expect(out.normal.approxEql(.init(1, 0, 0)));
    try testing.expectApproxEqAbs(@as(f32, -0.35), out.points[0].separation, 1e-5);
}

test "a cylinder standing on a box rests on its flat end, and two capsules side by side touch along their length" {
    var sa: Storage = .{};
    var sb: Storage = .{};
    const ground = solid(.{ .box = .{ .half = .init(5, 0.5, 5) } }, .at(.zero), &sa);
    const standing = solid(.{ .cylinder = .{ .half_height = 1, .radius = 0.5 } }, .at(.init(0, 1.49, 0)), &sb);
    const m = collide(ground, standing, 0.02);
    try testing.expectEqual(@as(u32, 4), m.count);
    try testing.expect(m.normal.approxEql(.init(0, 1, 0)));
    // The four spread over the round end.
    var spread: f32 = 0;
    for (m.pointSlice()) |p| spread = @max(spread, @abs(p.point.x) + @abs(p.point.z));
    try testing.expect(spread > 0.4);

    const a = solid(.{ .capsule = .{ .half_height = 1, .radius = 0.5 } }, .at(.zero), &sa);
    const b = solid(.{ .capsule = .{ .half_height = 1, .radius = 0.5 } }, .at(.init(0.95, 0.5, 0)), &sb);
    const side = collide(a, b, 0.02);
    try testing.expectEqual(@as(u32, 2), side.count);
    try testing.expect(side.normal.approxEql(.init(1, 0, 0)));
}

test "a box on a triangle of a floor touches it from above" {
    var sa: Storage = .{};
    var sb: Storage = .{};
    const tri = triangleSolid(&sa, .init(-5, 0, -5), .init(-5, 0, 5), .init(5, 0, 0));
    const crate = solid(.{ .box = .{ .half = .init(0.5, 0.5, 0.5) } }, .at(.init(0, 0.49, 0)), &sb);
    const m = collide(tri, crate, 0.02);
    try testing.expectEqual(@as(u32, 4), m.count);
    try testing.expect(m.normal.approxEql(.init(0, 1, 0)));
}

test "two hulls resting face to face touch at the corners of the smaller" {
    var points: [8]Vec3 = undefined;
    for (&points, 0..) |*v, i| v.* = .init(if (i & 1 != 0) 1 else -1, if (i & 2 != 0) 1 else -1, if (i & 4 != 0) 1 else -1);
    var cube = try hull_mod.init(testing.allocator, &points);
    defer cube.deinit(testing.allocator);
    var sa: Storage = .{};
    var sb: Storage = .{};
    const a = solid(.{ .hull = &cube }, .at(.zero), &sa);
    const b = solid(.{ .hull = &cube }, .{ .position = .init(0.5, 1.99, 0.3), .rotation = .fromAxisAngle(.init(0, 1, 0), 0.3) }, &sb);
    const m = collide(a, b, 0.02);
    try testing.expectEqual(@as(u32, 4), m.count);
    try testing.expect(m.normal.approxEql(.init(0, 1, 0)));
    for (m.pointSlice()) |p| try testing.expectApproxEqAbs(@as(f32, -0.01), p.separation, 1e-4);
}

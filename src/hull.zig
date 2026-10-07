// SPDX-License-Identifier: BSD-2-Clause

//! Convex solids made of flat faces: a hull built from any points, and the
//! box and the triangle seen the same way. What `collide` asks of a solid
//! with corners is a `Poly`: its corners, its faces - each a plane and its
//! corners round it, counter-clockwise seen from outside - and its edges.
//!
//! ```zig
//! var rock = try Hull.init(gpa, &points);
//! defer rock.deinit(gpa);
//! _ = try world.addShape(body, .{ .geometry = .{ .hull = &rock } });
//! ```
//!
//! **A hull is built one point at a time**: a first solid of four points,
//! and each point outside it after that seen from the faces it is outside
//! of, which go, and joined to the edge round them - the horizon - by new
//! triangles. Points inside are passed over. Then triangles that lie in one
//! plane are made one face: a box built from its eight corners has six
//! faces, not twelve triangles, which is what lets a box lying on another
//! touch it at four corners and rest flat.
//!
//! **Its mass is the sum of tetrahedra**, one from a point inside to each
//! triangle of each face, each with the inertia of a solid tetrahedron
//! (Blow and Binstock's covariance), moved to the hull's centre of mass.

const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Mat3 = geometry.Mat3;
const Aabb = geometry.Aabb;

/// One face: its plane, `normal . x = distance`, out of the solid, and
/// where its corners are in `face_vertices`.
pub const Face = struct {
    normal: Vec3,
    distance: f32,
    first: u16,
    count: u16,
};

/// A solid with corners, as `collide` sees it.
pub const Poly = struct {
    vertices: []const Vec3,
    faces: []const Face,
    face_vertices: []const u16,
    /// Each edge once - or, for a box and a cylinder, one edge each way
    /// they run: what the separating axes between two solids are crossed
    /// from.
    edges: []const [2]u16,
    /// For a hull, the two faces either side of each edge: what tells which
    /// pairs of edges can be nearest, so two hulls need not try every pair.
    /// Empty where `edges` stands for every edge running its way.
    edge_faces: []const [2]u16 = &.{},

    /// The corner furthest along `d`.
    pub fn support(self: Poly, d: Vec3) u16 {
        var best: u16 = 0;
        var most = self.vertices[0].dot(d);
        for (self.vertices[1..], 1..) |v, i| {
            const along = v.dot(d);
            if (along > most) {
                most = along;
                best = @intCast(i);
            }
        }
        return best;
    }

    pub fn cornersOf(self: Poly, face: Face) []const u16 {
        return self.face_vertices[face.first..][0..face.count];
    }
};

pub const Hull = @This();

vertices: []Vec3,
faces: []Face,
face_vertices: []u16,
edges: [][2]u16,
edge_faces: [][2]u16,
/// Where its mass is centred, and its volume.
centroid: Vec3,
volume: f32,
/// About the centroid, for a density of one.
inertia: Mat3,
bounds: Aabb,

pub const Error = Allocator.Error || error{
    /// Fewer than four points, or all of them in one plane or on one line:
    /// no solid.
    FlatHull,
};

pub fn poly(self: *const Hull) Poly {
    return .{ .vertices = self.vertices, .faces = self.faces, .face_vertices = self.face_vertices, .edges = self.edges, .edge_faces = self.edge_faces };
}

pub fn deinit(self: *Hull, gpa: Allocator) void {
    gpa.free(self.vertices);
    gpa.free(self.faces);
    gpa.free(self.face_vertices);
    gpa.free(self.edges);
    gpa.free(self.edge_faces);
    self.* = undefined;
}

/// The corner furthest along `d`, in the hull's own frame.
pub fn support(self: *const Hull, d: Vec3) Vec3 {
    return self.vertices[self.poly().support(d)];
}

/// The hull of `points`: every one of them on or inside it.
pub fn init(gpa: Allocator, points: []const Vec3) Error!Hull {
    if (points.len < 4) return error.FlatHull;
    const bounds = Aabb.ofPoints(points);
    const size = bounds.max.sub(bounds.min).maxComp();
    if (!(size > 0)) return error.FlatHull;
    // What counts as outside a face: a little past it, for the size.
    const eps = size * 1e-5;

    var builder: Builder = .{ .gpa = gpa, .points = points, .eps = eps };
    defer builder.deinit();
    try builder.start();
    for (points, 0..) |_, i| try builder.add(@intCast(i));
    return builder.finish(bounds);
}

/// A triangle of the hull as it is built: its corners, counter-clockwise
/// seen from outside, and its plane.
const Triangle = struct {
    v: [3]u32,
    normal: Vec3,
    distance: f32,
    alive: bool = true,
};

const Builder = struct {
    gpa: Allocator,
    points: []const Vec3,
    eps: f32,
    triangles: std.ArrayList(Triangle) = .empty,
    /// Each directed edge of a live triangle, to the triangle it is in.
    edge_owner: std.AutoHashMapUnmanaged([2]u32, u32) = .empty,
    used: std.ArrayList(bool) = .empty,

    fn deinit(self: *Builder) void {
        self.triangles.deinit(self.gpa);
        self.edge_owner.deinit(self.gpa);
        self.used.deinit(self.gpa);
    }

    fn makeTriangle(self: *Builder, a: u32, b: u32, c: u32) Allocator.Error!void {
        const pa = self.points[a];
        const n = self.points[b].sub(pa).cross(self.points[c].sub(pa));
        const unit = n.tryNorm() orelse n;
        const at: u32 = @intCast(self.triangles.items.len);
        try self.triangles.append(self.gpa, .{ .v = .{ a, b, c }, .normal = unit, .distance = unit.dot(pa) });
        try self.edge_owner.put(self.gpa, .{ a, b }, at);
        try self.edge_owner.put(self.gpa, .{ b, c }, at);
        try self.edge_owner.put(self.gpa, .{ c, a }, at);
    }

    fn kill(self: *Builder, at: u32) void {
        const t = &self.triangles.items[at];
        t.alive = false;
        for (0..3) |i| _ = self.edge_owner.remove(.{ t.v[i], t.v[(i + 1) % 3] });
    }

    /// A first solid of four points far apart: the two furthest along an
    /// axis, the furthest from the line between them, and the furthest from
    /// the plane of those three.
    fn start(self: *Builder) Error!void {
        try self.used.appendNTimes(self.gpa, false, self.points.len);
        const points = self.points;
        var lo: u32 = 0;
        var hi: u32 = 0;
        var best_span: f32 = -1;
        for ([_]Vec3{ .init(1, 0, 0), .init(0, 1, 0), .init(0, 0, 1) }) |axis| {
            var a: u32 = 0;
            var b: u32 = 0;
            for (points, 0..) |p, i| {
                if (p.dot(axis) < points[a].dot(axis)) a = @intCast(i);
                if (p.dot(axis) > points[b].dot(axis)) b = @intCast(i);
            }
            const span = points[b].dot(axis) - points[a].dot(axis);
            if (span > best_span) {
                best_span = span;
                lo = a;
                hi = b;
            }
        }
        const line = points[hi].sub(points[lo]);
        var third: u32 = lo;
        var furthest: f32 = 0;
        for (points, 0..) |p, i| {
            const away = line.cross(p.sub(points[lo])).lenSq();
            if (away > furthest) {
                furthest = away;
                third = @intCast(i);
            }
        }
        if (furthest <= self.eps * self.eps * line.lenSq()) return error.FlatHull;
        const plane = line.cross(points[third].sub(points[lo])).norm();
        var fourth: u32 = lo;
        var deepest: f32 = 0;
        for (points, 0..) |p, i| {
            const away = @abs(plane.dot(p.sub(points[lo])));
            if (away > deepest) {
                deepest = away;
                fourth = @intCast(i);
            }
        }
        if (deepest <= self.eps) return error.FlatHull;

        // Wound so every face looks out: the fourth point behind the first.
        var a = lo;
        var b = hi;
        const c = third;
        if (plane.dot(points[fourth].sub(points[lo])) > 0) std.mem.swap(u32, &a, &b);
        try self.makeTriangle(a, b, c);
        try self.makeTriangle(a, fourth, b);
        try self.makeTriangle(b, fourth, c);
        try self.makeTriangle(c, fourth, a);
        for ([_]u32{ a, b, c, fourth }) |i| self.used.items[i] = true;
    }

    /// Take `p` in, if it is outside.
    fn add(self: *Builder, p: u32) Allocator.Error!void {
        if (self.used.items[p]) return;
        const point = self.points[p];
        var visible: std.ArrayList(u32) = .empty;
        defer visible.deinit(self.gpa);
        for (self.triangles.items, 0..) |t, i| {
            if (t.alive and t.normal.dot(point) - t.distance > self.eps) try visible.append(self.gpa, @intCast(i));
        }
        if (visible.items.len == 0) return;
        self.used.items[p] = true;

        // The horizon: each edge of a face it sees whose other side it does
        // not see.
        var horizon: std.ArrayList([2]u32) = .empty;
        defer horizon.deinit(self.gpa);
        for (visible.items) |at| {
            const t = self.triangles.items[at];
            for (0..3) |i| {
                const a = t.v[i];
                const b = t.v[(i + 1) % 3];
                const across = self.edge_owner.get(.{ b, a }) orelse continue;
                if (std.mem.indexOfScalar(u32, visible.items, across) == null) try horizon.append(self.gpa, .{ a, b });
            }
        }
        for (visible.items) |at| self.kill(at);
        for (horizon.items) |edge| try self.makeTriangle(edge[0], edge[1], p);
    }

    /// The live triangles as faces - those in one plane merged - and the
    /// corners they use, renumbered.
    fn finish(self: *Builder, bounds: Aabb) Error!Hull {
        const gpa = self.gpa;
        var live: std.ArrayList(u32) = .empty;
        defer live.deinit(gpa);
        for (self.triangles.items, 0..) |t, i| if (t.alive) try live.append(gpa, @intCast(i));

        // Corners renumbered in the order they are first met.
        var renumber: std.AutoHashMapUnmanaged(u32, u16) = .empty;
        defer renumber.deinit(gpa);
        var vertices: std.ArrayList(Vec3) = .empty;
        errdefer vertices.deinit(gpa);
        for (live.items) |at| for (self.triangles.items[at].v) |v| {
            const entry = try renumber.getOrPut(gpa, v);
            if (!entry.found_existing) {
                entry.value_ptr.* = @intCast(vertices.items.len);
                try vertices.append(gpa, self.points[v]);
            }
        };

        // Faces: triangles flood-filled across edges while they stay in the
        // first one's plane.
        const grouped = try gpa.alloc(bool, self.triangles.items.len);
        defer gpa.free(grouped);
        @memset(grouped, false);
        var faces: std.ArrayList(Face) = .empty;
        errdefer faces.deinit(gpa);
        var face_vertices: std.ArrayList(u16) = .empty;
        errdefer face_vertices.deinit(gpa);
        var group: std.ArrayList(u32) = .empty;
        defer group.deinit(gpa);
        var boundary: std.AutoHashMapUnmanaged(u32, u32) = .empty;
        defer boundary.deinit(gpa);
        for (live.items) |seed| {
            if (grouped[seed]) continue;
            const plane = self.triangles.items[seed];
            group.clearRetainingCapacity();
            try group.append(gpa, seed);
            grouped[seed] = true;
            var next: usize = 0;
            while (next < group.items.len) : (next += 1) {
                const t = self.triangles.items[group.items[next]];
                for (0..3) |i| {
                    const across = self.edge_owner.get(.{ t.v[(i + 1) % 3], t.v[i] }) orelse continue;
                    if (grouped[across]) continue;
                    const other = self.triangles.items[across];
                    if (other.normal.dot(plane.normal) < 1 - 1e-4) continue;
                    const off = for (other.v) |v| {
                        if (@abs(plane.normal.dot(self.points[v]) - plane.distance) > self.eps * 4) break true;
                    } else false;
                    if (off) continue;
                    grouped[across] = true;
                    try group.append(gpa, across);
                }
            }
            // Its edge: each directed edge of the group whose reverse is not
            // in it, chained from one end to the next.
            boundary.clearRetainingCapacity();
            for (group.items) |at| {
                const t = self.triangles.items[at];
                for (0..3) |i| {
                    const a = t.v[i];
                    const b = t.v[(i + 1) % 3];
                    const across = self.edge_owner.get(.{ b, a }) orelse continue;
                    if (std.mem.indexOfScalar(u32, group.items, across) != null) continue;
                    try boundary.put(gpa, a, b);
                }
            }
            const first: u16 = @intCast(face_vertices.items.len);
            var it = boundary.iterator();
            const begin = (it.next() orelse continue).key_ptr.*;
            var at = begin;
            var count: u16 = 0;
            while (true) {
                try face_vertices.append(gpa, renumber.get(at).?);
                count += 1;
                at = boundary.get(at) orelse break;
                if (at == begin or count > boundary.count()) break;
            }
            // Averaged over its triangles, the plane is the face's own.
            var normal: Vec3 = .zero;
            for (group.items) |g| normal = normal.add(self.triangles.items[g].normal);
            normal = normal.norm();
            try faces.append(gpa, .{
                .normal = normal,
                .distance = normal.dot(vertices.items[face_vertices.items[first]]),
                .first = first,
                .count = count,
            });
        }

        // Each edge once, with the face either side of it.
        var edge_set: std.AutoArrayHashMapUnmanaged([2]u16, [2]u16) = .empty;
        defer edge_set.deinit(gpa);
        for (faces.items, 0..) |face, f| {
            const corners = face_vertices.items[face.first..][0..face.count];
            for (corners, 0..) |a, i| {
                const b = corners[(i + 1) % corners.len];
                const entry = try edge_set.getOrPut(gpa, if (a < b) .{ a, b } else .{ b, a });
                if (!entry.found_existing) {
                    entry.value_ptr.* = .{ @intCast(f), @intCast(f) };
                } else entry.value_ptr.*[1] = @intCast(f);
            }
        }
        const edges = try gpa.dupe([2]u16, edge_set.keys());
        errdefer gpa.free(edges);
        const edge_faces = try gpa.dupe([2]u16, edge_set.values());
        errdefer gpa.free(edge_faces);

        var made: Hull = .{
            .vertices = try vertices.toOwnedSlice(gpa),
            .faces = &.{},
            .face_vertices = &.{},
            .edges = edges,
            .edge_faces = edge_faces,
            .centroid = .zero,
            .volume = 0,
            .inertia = .zero,
            .bounds = bounds,
        };
        errdefer gpa.free(made.vertices);
        made.faces = try faces.toOwnedSlice(gpa);
        errdefer gpa.free(made.faces);
        made.face_vertices = try face_vertices.toOwnedSlice(gpa);
        made.bounds = Aabb.ofPoints(made.vertices);
        const mass = massOf(made.poly());
        made.centroid = mass.center;
        made.volume = mass.volume;
        made.inertia = mass.inertia;
        return made;
    }
};

/// What a solid of density one weighs, where, and how it turns.
pub const MassOf = struct {
    volume: f32,
    center: Vec3,
    /// About `center`.
    inertia: Mat3,
};

/// The mass of a closed solid with corners: a tetrahedron from a point
/// inside to each triangle of each face.
pub fn massOf(p: Poly) MassOf {
    var inside: Vec3 = .zero;
    for (p.vertices) |v| inside = inside.add(v);
    inside = inside.scale(1 / @as(f32, @floatFromInt(p.vertices.len)));

    var volume: f32 = 0;
    var moment: Vec3 = .zero;
    var covariance: Mat3 = .zero;
    for (p.faces) |face| {
        const corners = p.cornersOf(face);
        const a = p.vertices[corners[0]].sub(inside);
        for (1..corners.len - 1) |i| {
            const b = p.vertices[corners[i]].sub(inside);
            const c = p.vertices[corners[i + 1]].sub(inside);
            const det = a.dot(b.cross(c));
            const v = det / 6;
            volume += v;
            moment = moment.add(a.add(b).add(c).scale(v / 4));
            // A tetrahedron with one corner at the origin: its covariance.
            const sum = a.add(b).add(c);
            const cov = geometry.addMat(geometry.addMat(geometry.outer(a, a), geometry.outer(b, b)), geometry.addMat(geometry.outer(c, c), geometry.outer(sum, sum)));
            covariance = geometry.addMat(covariance, geometry.scaleMat(cov, det / 120));
        }
    }
    if (volume <= 0) return .{ .volume = 0, .center = inside, .inertia = .zero };
    const offset = moment.scale(1 / volume);
    // About the centre of mass, and from covariance to inertia.
    const about = geometry.addMat(covariance, geometry.scaleMat(geometry.outer(offset, offset), -volume));
    const trace = about.cols[0].x + about.cols[1].y + about.cols[2].z;
    const inertia = geometry.addMat(geometry.diagonal(.splat(trace)), geometry.scaleMat(about, -1));
    return .{ .volume = volume, .center = inside.add(offset), .inertia = inertia };
}

// -------------------------------------------------------------------------
// A box and a triangle, seen as solids with corners
// -------------------------------------------------------------------------

/// The corners are numbered by which way each is from the centre: bit 0 for
/// `+x`, bit 1 for `+y`, bit 2 for `+z`.
const box_faces = [_]Face{
    .{ .normal = .init(1, 0, 0), .distance = 0, .first = 0, .count = 4 },
    .{ .normal = .init(-1, 0, 0), .distance = 0, .first = 4, .count = 4 },
    .{ .normal = .init(0, 1, 0), .distance = 0, .first = 8, .count = 4 },
    .{ .normal = .init(0, -1, 0), .distance = 0, .first = 12, .count = 4 },
    .{ .normal = .init(0, 0, 1), .distance = 0, .first = 16, .count = 4 },
    .{ .normal = .init(0, 0, -1), .distance = 0, .first = 20, .count = 4 },
};
const box_face_vertices = [_]u16{ 1, 3, 7, 5, 0, 4, 6, 2, 2, 6, 7, 3, 0, 1, 5, 4, 4, 5, 7, 6, 0, 2, 3, 1 };
/// One edge each way: a box's edges run three ways only.
const box_edges = [_][2]u16{ .{ 0, 1 }, .{ 0, 2 }, .{ 0, 4 } };

/// A box of half size `half`, as a solid with corners, in its own frame.
pub const BoxPoly = struct {
    vertices: [8]Vec3,
    faces: [6]Face,

    pub fn init(half: Vec3) BoxPoly {
        var out: BoxPoly = .{ .vertices = undefined, .faces = box_faces };
        for (&out.vertices, 0..) |*v, i| v.* = .init(
            if (i & 1 != 0) half.x else -half.x,
            if (i & 2 != 0) half.y else -half.y,
            if (i & 4 != 0) half.z else -half.z,
        );
        for (&out.faces) |*f| f.distance = @abs(f.normal.dot(half));
        return out;
    }

    pub fn poly(self: *const BoxPoly) Poly {
        return .{ .vertices = &self.vertices, .faces = &self.faces, .face_vertices = &box_face_vertices, .edges = &box_edges };
    }
};

/// How many sides a cylinder has, colliding.
pub const cylinder_sides = 16;

/// A cylinder as a prism of `cylinder_sides` sides, its corners on the
/// round. Corners `0..16` go round the bottom and `16..32` the top.
pub const CylinderPoly = struct {
    vertices: [2 * cylinder_sides]Vec3,
    faces: [cylinder_sides + 2]Face,

    const face_vertices: [cylinder_sides * 6]u16 = blk: {
        var out: [cylinder_sides * 6]u16 = undefined;
        // The sides: up the one corner and down the next, out of the side.
        for (0..cylinder_sides) |i| {
            const next = (i + 1) % cylinder_sides;
            out[i * 4 + 0] = i;
            out[i * 4 + 1] = cylinder_sides + i;
            out[i * 4 + 2] = cylinder_sides + next;
            out[i * 4 + 3] = next;
        }
        // The bottom round the way that faces down, the top the other way.
        for (0..cylinder_sides) |i| {
            out[cylinder_sides * 4 + i] = i;
            out[cylinder_sides * 5 + i] = 2 * cylinder_sides - 1 - i;
        }
        break :blk out;
    };

    /// One edge each way: half the rim's edges, the other half running
    /// the same ways, and one up a side.
    const edges: [cylinder_sides / 2 + 1][2]u16 = blk: {
        var out: [cylinder_sides / 2 + 1][2]u16 = undefined;
        for (0..cylinder_sides / 2) |i| out[i] = .{ i, i + 1 };
        out[cylinder_sides / 2] = .{ 0, cylinder_sides };
        break :blk out;
    };

    pub fn init(half_height: f32, radius: f32) CylinderPoly {
        var out: CylinderPoly = .{ .vertices = undefined, .faces = undefined };
        const step = 2 * std.math.pi / @as(f32, cylinder_sides);
        for (0..cylinder_sides) |i| {
            const angle = step * @as(f32, @floatFromInt(i));
            const x = radius * @cos(angle);
            const z = radius * @sin(angle);
            out.vertices[i] = .init(x, -half_height, z);
            out.vertices[cylinder_sides + i] = .init(x, half_height, z);
            const middle = angle + step / 2;
            const n = Vec3.init(@cos(middle), 0, @sin(middle));
            out.faces[i] = .{ .normal = n, .distance = radius * @cos(step / 2), .first = @intCast(i * 4), .count = 4 };
        }
        out.faces[cylinder_sides] = .{ .normal = .init(0, -1, 0), .distance = half_height, .first = cylinder_sides * 4, .count = cylinder_sides };
        out.faces[cylinder_sides + 1] = .{ .normal = .init(0, 1, 0), .distance = half_height, .first = cylinder_sides * 5, .count = cylinder_sides };
        return out;
    }

    pub fn poly(self: *const CylinderPoly) Poly {
        return .{ .vertices = &self.vertices, .faces = &self.faces, .face_vertices = &face_vertices, .edges = &edges };
    }
};

const triangle_face_vertices = [_]u16{ 0, 1, 2, 0, 2, 1 };
const triangle_edges = [_][2]u16{ .{ 0, 1 }, .{ 1, 2 }, .{ 2, 0 } };

/// A triangle, as a solid with no thickness: a face each way.
pub const TrianglePoly = struct {
    vertices: [3]Vec3,
    faces: [2]Face,

    pub fn init(a: Vec3, b: Vec3, c: Vec3) TrianglePoly {
        const n = b.sub(a).cross(c.sub(a)).tryNorm() orelse Vec3.init(0, 1, 0);
        return .{
            .vertices = .{ a, b, c },
            .faces = .{
                .{ .normal = n, .distance = n.dot(a), .first = 0, .count = 3 },
                .{ .normal = n.neg(), .distance = -n.dot(a), .first = 3, .count = 3 },
            },
        };
    }

    pub fn poly(self: *const TrianglePoly) Poly {
        return .{ .vertices = &self.vertices, .faces = &self.faces, .face_vertices = &triangle_face_vertices, .edges = &triangle_edges };
    }
};

// -------------------------------------------------------------------------
// Tests
// -------------------------------------------------------------------------

fn cubeCorners() [8]Vec3 {
    var out: [8]Vec3 = undefined;
    for (&out, 0..) |*v, i| v.* = .init(if (i & 1 != 0) 1 else -1, if (i & 2 != 0) 1 else -1, if (i & 4 != 0) 1 else -1);
    return out;
}

test "a cube's corners make six square faces, twelve edges, and a cube's mass" {
    var points: [20]Vec3 = undefined;
    @memcpy(points[0..8], &cubeCorners());
    // Points inside and on faces, which the hull passes over.
    for (points[8..], 0..) |*p, i| p.* = .init(0.1 * @as(f32, @floatFromInt(i)) - 0.6, 0.05, if (i % 2 == 0) 1 else 0);
    var cube = try Hull.init(testing.allocator, &points);
    defer cube.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 8), cube.vertices.len);
    try testing.expectEqual(@as(usize, 6), cube.faces.len);
    try testing.expectEqual(@as(usize, 12), cube.edges.len);
    // Each edge between two faces at right angles.
    for (cube.edge_faces) |faces| try testing.expectApproxEqAbs(@as(f32, 0), cube.faces[faces[0]].normal.dot(cube.faces[faces[1]].normal), 1e-5);
    for (cube.faces) |f| {
        try testing.expectEqual(@as(u16, 4), f.count);
        try testing.expectApproxEqAbs(@as(f32, 1), f.distance, 1e-5);
        // Every corner of the face on its plane.
        for (cube.poly().cornersOf(f)) |v| try testing.expectApproxEqAbs(f.distance, f.normal.dot(cube.vertices[v]), 1e-5);
    }
    try testing.expectApproxEqAbs(@as(f32, 8), cube.volume, 1e-4);
    try testing.expect(cube.centroid.approxEql(.zero));
    // A solid cube of side two: m (a^2 + a^2) / 12 = 8 * 8 / 12.
    try testing.expectApproxEqAbs(@as(f32, 64.0 / 12.0), cube.inertia.cols[0].x, 1e-3);
    try testing.expectApproxEqAbs(@as(f32, 0), cube.inertia.cols[0].y, 1e-4);
}

test "the mass of a box seen as a solid with corners is the box's" {
    const box = BoxPoly.init(.init(1, 2, 3));
    const mass = massOf(box.poly());
    try testing.expectApproxEqAbs(@as(f32, 48), mass.volume, 1e-3);
    // Ixx = m (b^2 + c^2) / 3 for half sizes a, b, c.
    try testing.expectApproxEqAbs(@as(f32, 48 * (4 + 9) / 3.0), mass.inertia.cols[0].x, 1e-2);
    try testing.expectApproxEqAbs(@as(f32, 48 * (1 + 4) / 3.0), mass.inertia.cols[2].z, 1e-2);
    // Each face's corners wound out of it.
    const p = box.poly();
    for (p.faces) |f| {
        const c = p.cornersOf(f);
        const n = p.vertices[c[1]].sub(p.vertices[c[0]]).cross(p.vertices[c[2]].sub(p.vertices[c[1]])).norm();
        try testing.expect(n.approxEql(f.normal));
    }
}

test "a cylinder's prism has its faces wound out of it and its corners on its round" {
    const c = CylinderPoly.init(2, 1);
    const p = c.poly();
    for (p.faces) |f| {
        const corners = p.cornersOf(f);
        const n = p.vertices[corners[1]].sub(p.vertices[corners[0]]).cross(p.vertices[corners[2]].sub(p.vertices[corners[1]])).norm();
        try testing.expect(n.approxEql(f.normal));
        for (corners) |v| try testing.expectApproxEqAbs(f.distance, f.normal.dot(p.vertices[v]), 1e-5);
    }
    const mass = massOf(p);
    // A little less than the cylinder round it.
    try testing.expect(mass.volume < 4 * std.math.pi and mass.volume > 0.97 * 4 * std.math.pi);
}

test "points all in one plane make no solid, and a hull of a sphere's points holds them all" {
    const flat = [_]Vec3{ .init(0, 0, 0), .init(1, 0, 0), .init(0, 0, 1), .init(1, 0, 1), .init(0.5, 0, 0.5) };
    try testing.expectError(error.FlatHull, Hull.init(testing.allocator, &flat));

    var points: [200]Vec3 = undefined;
    var random = std.Random.DefaultPrng.init(7);
    for (&points) |*p| p.* = Vec3.init(random.random().floatNorm(f32), random.random().floatNorm(f32), random.random().floatNorm(f32)).norm();
    var ball = try Hull.init(testing.allocator, &points);
    defer ball.deinit(testing.allocator);
    for (points) |p| for (ball.faces) |f| try testing.expect(f.normal.dot(p) - f.distance < 1e-4);
    // A sphere of radius one is 4.19; its hull is a little less.
    try testing.expect(ball.volume > 3.6 and ball.volume < 4.19);
}

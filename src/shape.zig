// SPDX-License-Identifier: BSD-2-Clause

//! What a body collides with: a ball, a box, a capsule, a cylinder, a
//! convex hull or a mesh of triangles, placed on the body, with what it is
//! made of and what it touches.
//!
//! ```zig
//! _ = try world.addShape(crate, .box(.init(0.5, 0.5, 0.5)));
//! _ = try world.addShape(player, .{ .geometry = .{ .capsule = .{ .half_height = 0.6, .radius = 0.3 } }, .material = .{ .friction = 0 } });
//! ```
//!
//! **A capsule and a cylinder stand along their own `y`**, as a character
//! stands; turn the shape on its body to lay one down. A capsule is every
//! point within its radius of a segment, so it slides over a step's edge;
//! a cylinder is flat at both ends and, colliding, a prism of sixteen
//! sides - round enough to roll, flat enough to stand.
//!
//! **Mass comes from density and volume**, so a big crate weighs more than a
//! small one without anyone typing a mass. A mesh has no inside and no
//! mass: it is for a static body.

const std = @import("std");
const testing = std.testing;

const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Mat3 = geometry.Mat3;
const Transform = geometry.Transform;
const Aabb = geometry.Aabb;
const Hull = @import("hull.zig");
const Mesh = @import("mesh.zig");

pub const Sphere = struct {
    radius: f32,
};

pub const Box = struct {
    /// Half its size each way.
    half: Vec3,
};

/// Every point within `radius` of the segment from `-half_height` to
/// `half_height` along `y`.
pub const Capsule = struct {
    half_height: f32,
    radius: f32,
};

/// Round about `y`, flat at `-half_height` and `half_height`.
pub const Cylinder = struct {
    half_height: f32,
    radius: f32,
};

pub const Geometry = union(enum) {
    sphere: Sphere,
    box: Box,
    capsule: Capsule,
    cylinder: Cylinder,
    /// Kept by the caller for as long as a shape uses it.
    hull: *const Hull,
    /// Kept by the caller for as long as a shape uses it; on a static body.
    mesh: *const Mesh,

    /// The box round it, placed by `xf`.
    pub fn aabb(self: Geometry, xf: Transform) Aabb {
        return switch (self) {
            .sphere => |s| .{ .min = xf.position.sub(.splat(s.radius)), .max = xf.position.add(.splat(s.radius)) },
            .capsule => |c| blk: {
                const axis = xf.turn(.init(0, c.half_height, 0));
                const a = xf.position.sub(axis);
                const b = xf.position.add(axis);
                break :blk (Aabb{ .min = a.min(b), .max = a.max(b) }).grow(c.radius);
            },
            .box => |b| turnedBox(xf, .{ .min = b.half.neg(), .max = b.half }),
            .cylinder => |c| turnedBox(xf, .{ .min = .init(-c.radius, -c.half_height, -c.radius), .max = .init(c.radius, c.half_height, c.radius) }),
            .hull => |h| turnedBox(xf, h.bounds),
            .mesh => |m| turnedBox(xf, m.bounds),
        };
    }

    /// What it weighs at `density`, where, and how it turns, in its own
    /// frame.
    pub fn massOf(self: Geometry, density: f32) Mass {
        const pi = std.math.pi;
        return switch (self) {
            .sphere => |s| blk: {
                const m = density * 4.0 / 3.0 * pi * s.radius * s.radius * s.radius;
                break :blk .{ .mass = m, .center = .zero, .inertia = geometry.diagonal(.splat(0.4 * m * s.radius * s.radius)) };
            },
            .box => |b| blk: {
                const m = density * 8 * b.half.x * b.half.y * b.half.z;
                const x2 = b.half.x * b.half.x;
                const y2 = b.half.y * b.half.y;
                const z2 = b.half.z * b.half.z;
                break :blk .{ .mass = m, .center = .zero, .inertia = geometry.diagonal(Vec3.init(y2 + z2, x2 + z2, x2 + y2).scale(m / 3)) };
            },
            .capsule => |c| blk: {
                const r = c.radius;
                const h = c.half_height;
                const tube = density * pi * r * r * 2 * h;
                const ball = density * 4.0 / 3.0 * pi * r * r * r;
                const along = tube * r * r / 2 + ball * 0.4 * r * r;
                const across = tube * (h * h / 3 + r * r / 4) + ball * (0.4 * r * r + h * h + 0.75 * h * r);
                break :blk .{ .mass = tube + ball, .center = .zero, .inertia = geometry.diagonal(.init(across, along, across)) };
            },
            .cylinder => |c| blk: {
                const m = density * pi * c.radius * c.radius * 2 * c.half_height;
                const across = m * (c.radius * c.radius / 4 + c.half_height * c.half_height / 3);
                break :blk .{ .mass = m, .center = .zero, .inertia = geometry.diagonal(.init(across, m * c.radius * c.radius / 2, across)) };
            },
            .hull => |h| .{ .mass = density * h.volume, .center = h.centroid, .inertia = geometry.scaleMat(h.inertia, density) },
            .mesh => .{ .mass = 0, .center = .zero, .inertia = .zero },
        };
    }

    pub const RayHit = struct {
        fraction: f32,
        /// Out of the surface, back along the ray.
        normal: Vec3,
    };

    /// Where a ray from `origin` along `translation`, outside, first meets
    /// it - no further than `max_fraction` - placed by `xf`. A ray starting
    /// inside a solid meets nothing.
    pub fn rayCast(self: Geometry, xf: Transform, origin: Vec3, translation: Vec3, max_fraction: f32) ?RayHit {
        const o = xf.unapply(origin);
        const d = xf.unturn(translation);
        const local: ?RayHit = switch (self) {
            .sphere => |s| raySphere(o, d, .zero, s.radius),
            .box => |b| rayBox(o, d, b.half),
            .capsule => |c| rayCapsule(o, d, c.half_height, c.radius),
            .cylinder => |c| rayCylinder(o, d, c.half_height, c.radius),
            .hull => |h| rayHull(h, o, d),
            .mesh => |m| if (m.rayCast(o, d, max_fraction)) |hit| .{ .fraction = hit.fraction, .normal = hit.normal } else null,
        };
        const hit = local orelse return null;
        if (hit.fraction > max_fraction) return null;
        return .{ .fraction = hit.fraction, .normal = xf.turn(hit.normal) };
    }
};

/// The box round a box of the frame `xf` places, outside it.
fn turnedBox(xf: Transform, box: Aabb) Aabb {
    const r = geometry.rotationOf(xf.rotation);
    const center = xf.apply(box.center());
    const e = box.extent();
    const reach = Vec3.init(
        @abs(r.cols[0].x) * e.x + @abs(r.cols[1].x) * e.y + @abs(r.cols[2].x) * e.z,
        @abs(r.cols[0].y) * e.x + @abs(r.cols[1].y) * e.y + @abs(r.cols[2].y) * e.z,
        @abs(r.cols[0].z) * e.x + @abs(r.cols[1].z) * e.y + @abs(r.cols[2].z) * e.z,
    );
    return .{ .min = center.sub(reach), .max = center.add(reach) };
}

/// What a geometry weighs, in its own frame.
pub const Mass = struct {
    mass: f32,
    center: Vec3,
    /// About `center`.
    inertia: Mat3,
};

/// How bodies made of it rub and bounce, and how heavy it is.
pub const Material = struct {
    /// Nought is ice; one is rubber on concrete.
    friction: f32 = 0.6,
    /// How much of the speed a hit gives back: nought a beanbag, one a
    /// superball.
    restitution: f32 = 0,
    /// Mass a cubic unit.
    density: f32 = 1,
};

/// What it is, and what it looks for: two shapes touch when the rule
/// says the layers and the masks agree.
pub const Filter = struct {
    layer: u32 = 1,
    mask: u32 = 0xFFFF_FFFF,
};

/// How two filters agree.
pub const FilterRule = enum {
    /// Each one's mask has the other's layer.
    both,
    /// Either one's mask has the other's layer.
    either,

    pub fn touches(rule: FilterRule, a: Filter, b: Filter) bool {
        const ab = a.mask & b.layer != 0;
        const ba = b.mask & a.layer != 0;
        return switch (rule) {
            .both => ab and ba,
            .either => ab or ba,
        };
    }
};

/// How two surfaces' numbers make the pair's.
pub const Mix = enum {
    geometric_mean,
    average,
    minimum,
    maximum,
    multiply,
    /// The two added, no more than one: a bounce off something bouncy.
    sum_clamped,

    pub fn of(mix: Mix, a: f32, b: f32) f32 {
        return switch (mix) {
            .geometric_mean => @sqrt(@max(a * b, 0)),
            .average => (a + b) / 2,
            .minimum => @min(a, b),
            .maximum => @max(a, b),
            .multiply => a * b,
            .sum_clamped => @min(a + b, 1),
        };
    }
};

/// What `World.addShape` takes.
pub const Shape = struct {
    geometry: Geometry,
    /// Where it sits on its body, from the body's origin.
    offset: Transform = .identity,
    material: Material = .{},
    filter: Filter = .{},
    /// Says what overlaps it and pushes nothing: a trigger, a pickup.
    sensor: bool = false,
    /// The caller's: an entity, an index.
    user_data: u64 = 0,

    pub fn sphere(radius: f32) Shape {
        return .{ .geometry = .{ .sphere = .{ .radius = radius } } };
    }

    pub fn box(half: Vec3) Shape {
        return .{ .geometry = .{ .box = .{ .half = half } } };
    }

    pub fn capsule(half_height: f32, radius: f32) Shape {
        return .{ .geometry = .{ .capsule = .{ .half_height = half_height, .radius = radius } } };
    }

    pub fn cylinder(half_height: f32, radius: f32) Shape {
        return .{ .geometry = .{ .cylinder = .{ .half_height = half_height, .radius = radius } } };
    }
};

// -------------------------------------------------------------------------
// Rays, each in the shape's own frame
// -------------------------------------------------------------------------

fn raySphere(o: Vec3, d: Vec3, center: Vec3, radius: f32) ?Geometry.RayHit {
    const m = o.sub(center);
    const a = d.dot(d);
    if (a < 1e-20) return null;
    const b = m.dot(d);
    const c = m.dot(m) - radius * radius;
    // Inside, or past it going away.
    if (c <= 0) return null;
    if (b > 0) return null;
    const disc = b * b - a * c;
    if (disc < 0) return null;
    const t = (-b - @sqrt(disc)) / a;
    if (t < 0) return null;
    return .{ .fraction = t, .normal = m.mulAdd(d, t).norm() };
}

fn rayBox(o: Vec3, d: Vec3, half: Vec3) ?Geometry.RayHit {
    const lo = [3]f32{ -half.x, -half.y, -half.z };
    const hi = [3]f32{ half.x, half.y, half.z };
    const oo = [3]f32{ o.x, o.y, o.z };
    const dd = [3]f32{ d.x, d.y, d.z };
    var enter: f32 = -std.math.inf(f32);
    var leave: f32 = std.math.inf(f32);
    var axis: usize = 0;
    var sign: f32 = 1;
    for (0..3) |i| {
        if (@abs(dd[i]) < 1e-12) {
            if (oo[i] < lo[i] or oo[i] > hi[i]) return null;
            continue;
        }
        const inv = 1 / dd[i];
        var t1 = (lo[i] - oo[i]) * inv;
        var t2 = (hi[i] - oo[i]) * inv;
        var s: f32 = -1;
        if (t1 > t2) {
            std.mem.swap(f32, &t1, &t2);
            s = 1;
        }
        if (t1 > enter) {
            enter = t1;
            axis = i;
            sign = s;
        }
        leave = @min(leave, t2);
        if (enter > leave) return null;
    }
    if (enter < 0) return null;
    var n: [3]f32 = .{ 0, 0, 0 };
    n[axis] = sign;
    return .{ .fraction = enter, .normal = .init(n[0], n[1], n[2]) };
}

fn rayCapsule(o: Vec3, d: Vec3, half_height: f32, radius: f32) ?Geometry.RayHit {
    var best: ?Geometry.RayHit = null;
    // The side: a tube about y, between the two ends.
    const a = d.x * d.x + d.z * d.z;
    const b = o.x * d.x + o.z * d.z;
    const c = o.x * o.x + o.z * o.z - radius * radius;
    if (a > 1e-20 and c > 0 and b < 0) {
        const disc = b * b - a * c;
        if (disc >= 0) {
            const t = (-b - @sqrt(disc)) / a;
            const y = o.y + t * d.y;
            if (t >= 0 and y >= -half_height and y <= half_height) {
                best = .{ .fraction = t, .normal = Vec3.init(o.x + t * d.x, 0, o.z + t * d.z).norm() };
            }
        }
    }
    for ([_]f32{ -half_height, half_height }) |y| {
        const hit = raySphere(o, d, .init(0, y, 0), radius) orelse continue;
        if (best == null or hit.fraction < best.?.fraction) best = hit;
    }
    // Inside it is nothing to hit; a ray that starts in the tube between the
    // balls met their back sides only if it started outside them.
    const nearest = std.math.clamp(o.y, -half_height, half_height);
    if (o.sub(.init(0, nearest, 0)).lenSq() <= radius * radius) return null;
    return best;
}

fn rayCylinder(o: Vec3, d: Vec3, half_height: f32, radius: f32) ?Geometry.RayHit {
    if (@abs(o.y) <= half_height and o.x * o.x + o.z * o.z <= radius * radius) return null;
    var best: ?Geometry.RayHit = null;
    const a = d.x * d.x + d.z * d.z;
    const b = o.x * d.x + o.z * d.z;
    const c = o.x * o.x + o.z * o.z - radius * radius;
    if (a > 1e-20 and c > 0 and b < 0) {
        const disc = b * b - a * c;
        if (disc >= 0) {
            const t = (-b - @sqrt(disc)) / a;
            const y = o.y + t * d.y;
            if (t >= 0 and @abs(y) <= half_height) best = .{ .fraction = t, .normal = Vec3.init(o.x + t * d.x, 0, o.z + t * d.z).norm() };
        }
    }
    // The flat ends.
    if (@abs(d.y) > 1e-12) for ([_]f32{ -half_height, half_height }) |y| {
        const t = (y - o.y) / d.y;
        if (t < 0) continue;
        // Only from outside the end.
        if ((y > 0 and o.y < y) or (y < 0 and o.y > y)) continue;
        const x = o.x + t * d.x;
        const z = o.z + t * d.z;
        if (x * x + z * z > radius * radius) continue;
        if (best == null or t < best.?.fraction) best = .{ .fraction = t, .normal = .init(0, if (y > 0) 1 else -1, 0) };
    };
    return best;
}

/// A ray against a hull's planes: in through the last face it enters, out
/// through the first it leaves.
fn rayHull(h: *const Hull, o: Vec3, d: Vec3) ?Geometry.RayHit {
    var enter: f32 = 0;
    var leave: f32 = std.math.inf(f32);
    var normal: ?Vec3 = null;
    for (h.faces) |f| {
        const num = f.distance - f.normal.dot(o);
        const den = f.normal.dot(d);
        if (@abs(den) < 1e-12) {
            if (num < 0) return null;
            continue;
        }
        const t = num / den;
        if (den < 0) {
            if (t > enter) {
                enter = t;
                normal = f.normal;
            }
        } else leave = @min(leave, t);
        if (enter > leave) return null;
    }
    // A ray that starts inside never crosses a face going in.
    return .{ .fraction = enter, .normal = normal orelse return null };
}

test "a shape's box holds it however it is turned" {
    const turned: Transform = .{ .position = .init(1, 2, 3), .rotation = .fromAxisAngle(.init(0, 0, 1), std.math.pi / 4.0) };
    const box = (Geometry{ .box = .{ .half = .init(1, 1, 1) } }).aabb(turned);
    try testing.expectApproxEqAbs(@as(f32, 1 - std.math.sqrt2), box.min.x, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, 2), box.max.z - 2, 1e-5);
    const lying = (Geometry{ .capsule = .{ .half_height = 1, .radius = 0.5 } }).aabb(.{ .rotation = .fromAxisAngle(.init(0, 0, 1), std.math.pi / 2.0) });
    try testing.expectApproxEqAbs(@as(f32, 1.5), lying.max.x, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, 0.5), lying.max.y, 1e-5);
}

test "each shape's mass is its density times its volume" {
    const ball = (Geometry{ .sphere = .{ .radius = 1 } }).massOf(2);
    try testing.expectApproxEqAbs(@as(f32, 8.0 / 3.0 * std.math.pi), ball.mass, 1e-4);
    const capsule = (Geometry{ .capsule = .{ .half_height = 1, .radius = 1 } }).massOf(1);
    // A tube two long and a whole ball.
    try testing.expectApproxEqAbs(@as(f32, 2 * std.math.pi + 4.0 / 3.0 * std.math.pi), capsule.mass, 1e-4);
    // Harder to turn about a lying axis than about its own.
    try testing.expect(capsule.inertia.cols[0].x > capsule.inertia.cols[1].y);
    const cylinder = (Geometry{ .cylinder = .{ .half_height = 1, .radius = 1 } }).massOf(1);
    try testing.expectApproxEqAbs(@as(f32, std.math.pi), cylinder.inertia.cols[1].y, 1e-4);
}

test "rays meet each shape where its surface is, with the normal out of it" {
    const along = Vec3.init(10, 0, 0);
    const from = Vec3.init(-5, 0, 0);
    const placed: Transform = .at(.init(0, 0, 0));
    const ball = (Geometry{ .sphere = .{ .radius = 1 } }).rayCast(placed, from, along, 1).?;
    try testing.expectApproxEqAbs(@as(f32, 0.4), ball.fraction, 1e-6);
    try testing.expect(ball.normal.approxEql(.init(-1, 0, 0)));
    const box = (Geometry{ .box = .{ .half = .init(2, 1, 1) } }).rayCast(placed, from, along, 1).?;
    try testing.expectApproxEqAbs(@as(f32, 0.3), box.fraction, 1e-6);
    // A capsule lying along x: hit on its round end.
    const lying: Transform = .{ .rotation = .fromAxisAngle(.init(0, 0, 1), std.math.pi / 2.0) };
    const capsule = (Geometry{ .capsule = .{ .half_height = 1, .radius = 1 } }).rayCast(lying, from, along, 1).?;
    try testing.expectApproxEqAbs(@as(f32, 0.3), capsule.fraction, 1e-5);
    // Down onto a standing cylinder's top.
    const top = (Geometry{ .cylinder = .{ .half_height = 1, .radius = 1 } }).rayCast(placed, .init(0.5, 5, 0), .init(0, -10, 0), 1).?;
    try testing.expectApproxEqAbs(@as(f32, 0.4), top.fraction, 1e-6);
    try testing.expect(top.normal.approxEql(.init(0, 1, 0)));
    // From inside: nothing.
    try testing.expect((Geometry{ .sphere = .{ .radius = 1 } }).rayCast(placed, .zero, along, 1) == null);
}

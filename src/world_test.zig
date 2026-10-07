// SPDX-License-Identifier: BSD-2-Clause

//! Whole worlds, stepped: things falling, resting, stacking, rolling,
//! sleeping, and asked about.

const std = @import("std");
const testing = std.testing;

const World = @import("World.zig");
const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Quat = geometry.Quat;
const Hull = @import("hull.zig");
const Mesh = @import("mesh.zig");

const dt = 1.0 / 60.0;

fn stepFor(world: *World, seconds: f32) !void {
    var t: f32 = 0;
    while (t < seconds) : (t += dt) try world.step(dt);
}

/// A world with a floor whose top is at nought.
fn withFloor() !World {
    var world: World = .init(testing.allocator, .{});
    errdefer world.deinit();
    const ground = try world.createBody(.{ .type = .static, .position = .init(0, -0.5, 0) });
    _ = try world.addShape(ground, .box(.init(20, 0.5, 20)));
    return world;
}

test "a box dropped on the floor comes to rest on it, and falls asleep" {
    var world = try withFloor();
    defer world.deinit();
    const crate = try world.createBody(.{ .position = .init(0, 3, 0), .rotation = .fromAxisAngle(.init(0.3, 1, 0.2), 0.4) });
    _ = try world.addShape(crate, .box(.init(0.5, 0.5, 0.5)));
    try stepFor(&world, 4);
    const b = world.body(crate).?;
    // Flat on a face: its centre half a metre up, sunk no more than the slop.
    try testing.expectApproxEqAbs(@as(f32, 0.5), b.center.y, 0.02);
    try testing.expect(b.linear_velocity.len() < 0.05);
    try testing.expect(!b.awake);
    try testing.expectEqual(@as(usize, 0), world.awakeCount());
}

test "a stack of six boxes stands" {
    var world = try withFloor();
    defer world.deinit();
    var crates: [6]World.BodyId = undefined;
    for (&crates, 0..) |*c, i| {
        c.* = try world.createBody(.{ .position = .init(0, 0.5 + @as(f32, @floatFromInt(i)) * 1.0, 0) });
        _ = try world.addShape(c.*, .box(.init(0.5, 0.5, 0.5)));
    }
    try stepFor(&world, 5);
    const top = world.body(crates[5]).?;
    try testing.expectApproxEqAbs(@as(f32, 5.5), top.center.y, 0.1);
    try testing.expect(@abs(top.center.x) < 0.05 and @abs(top.center.z) < 0.05);
}

test "a ball rolls down a slope, turning as it goes" {
    var world = try withFloor();
    defer world.deinit();
    const slope = try world.createBody(.{ .type = .static, .position = .init(0, 2, 0), .rotation = .fromAxisAngle(.init(0, 0, 1), -0.3) });
    _ = try world.addShape(slope, .box(.init(5, 0.2, 2)));
    const ball = try world.createBody(.{ .position = .init(-3, 4, 0) });
    _ = try world.addShape(ball, .sphere(0.4));
    try stepFor(&world, 1.5);
    const b = world.body(ball).?;
    // Down the slope, towards +x, and spinning about -z as it rolls.
    try testing.expect(b.linear_velocity.x > 1);
    try testing.expect(b.angular_velocity.z < -1);
}

test "a box rests on a floor of triangles, and a capsule that cannot turn stands on it" {
    var world: World = .init(testing.allocator, .{});
    defer world.deinit();
    const vertices = [_]Vec3{ .init(-10, 0, -10), .init(10, 0, -10), .init(10, 0, 10), .init(-10, 0, 10) };
    var floor = try Mesh.init(testing.allocator, &vertices, &.{ 0, 2, 1, 0, 3, 2 });
    defer floor.deinit(testing.allocator);
    const ground = try world.createBody(.{ .type = .static });
    _ = try world.addShape(ground, .{ .geometry = .{ .mesh = &floor } });
    const crate = try world.createBody(.{ .position = .init(1, 2, 1) });
    _ = try world.addShape(crate, .box(.init(0.5, 0.5, 0.5)));
    const player = try world.createBody(.{ .position = .init(-2, 2, 0), .lock_rotation = true });
    _ = try world.addShape(player, .capsule(0.5, 0.3));
    try stepFor(&world, 3);
    try testing.expectApproxEqAbs(@as(f32, 0.5), world.body(crate).?.center.y, 0.02);
    try testing.expectApproxEqAbs(@as(f32, 0.8), world.body(player).?.center.y, 0.02);
    // A mesh on a body that moves is refused.
    const moving = try world.createBody(.{});
    try testing.expectError(error.MeshOnMovingBody, world.addShape(moving, .{ .geometry = .{ .mesh = &floor } }));
}

test "a ball falling through a sensor is seen going in and coming out, and pushed by nothing" {
    var world = try withFloor();
    defer world.deinit();
    const area = try world.createBody(.{ .type = .static, .position = .init(0, 3, 0) });
    const sensor = try world.addShape(area, .{ .geometry = .{ .box = .{ .half = .init(1, 0.5, 1) } }, .sensor = true });
    const ball = try world.createBody(.{ .position = .init(0, 6, 0) });
    _ = try world.addShape(ball, .sphere(0.25));
    var began = false;
    var ended = false;
    var t: f32 = 0;
    while (t < 2) : (t += dt) {
        try world.step(dt);
        for (world.beginEvents()) |e| if (e.sensor and (e.shape_a.eql(sensor) or e.shape_b.eql(sensor))) {
            began = true;
        };
        for (world.endEvents()) |e| if (e.sensor and (e.shape_a.eql(sensor) or e.shape_b.eql(sensor))) {
            ended = true;
        };
    }
    try testing.expect(began and ended);
    try testing.expect(world.body(ball).?.center.y < 0.3);
}

test "a ray meets the first thing in its way, and a cast stops where a shape first touches" {
    var world = try withFloor();
    defer world.deinit();
    const crate = try world.createBody(.{ .type = .static, .position = .init(0, 1, 0) });
    _ = try world.addShape(crate, .box(.init(1, 1, 1)));
    const hit = world.castRay(.init(0, 10, 0), .init(0, -20, 0), .{}).?;
    try testing.expectApproxEqAbs(@as(f32, 8.0 / 20.0), hit.fraction, 1e-5);
    try testing.expect(hit.body.eql(crate));
    try testing.expect(hit.normal.approxEql(.init(0, 1, 0)));
    // Passing over the crate: the floor.
    const past = world.castRay(.init(0, 10, 0), .init(0, -20, 0), .{ .ignore = crate }).?;
    try testing.expectApproxEqAbs(@as(f32, 0.5), past.fraction, 1e-5);

    // A ball cast down onto the crate stops its radius above it.
    const cast = world.castShape(.{ .sphere = .{ .radius = 0.5 } }, .at(.init(0.2, 6, 0.1)), .init(0, -10, 0), .{}).?;
    try testing.expectApproxEqAbs(@as(f32, (6 - 2.5) / 10.0), cast.fraction, 2e-3);
    try testing.expect(cast.normal.approxEql(.init(0, 1, 0)));
    // Sunk in: how deep, and which way out.
    var found: [4]World.Penetration = undefined;
    const deep = world.penetrations(.{ .box = .{ .half = .init(0.5, 0.5, 0.5) } }, .at(.init(0, 2.3, 0)), .{}, &found);
    try testing.expectEqual(@as(usize, 1), deep.len);
    try testing.expectApproxEqAbs(@as(f32, 0.2), deep[0].depth, 1e-4);
    try testing.expect(deep[0].normal.approxEql(.init(0, 1, 0)));
}

test "a box on a rising platform rides it up, and a bouncy ball bounces" {
    var world = try withFloor();
    defer world.deinit();
    const lift = try world.createBody(.{ .type = .kinematic, .position = .init(0, 1, 0), .linear_velocity = .init(0, 1, 0) });
    _ = try world.addShape(lift, .box(.init(1, 0.1, 1)));
    const crate = try world.createBody(.{ .position = .init(0, 1.6, 0) });
    _ = try world.addShape(crate, .box(.init(0.25, 0.25, 0.25)));
    try stepFor(&world, 2);
    // The lift is at three: the crate on it.
    try testing.expectApproxEqAbs(@as(f32, 3.35), world.body(crate).?.center.y, 0.05);

    const ball = try world.createBody(.{ .position = .init(5, 3, 5) });
    _ = try world.addShape(ball, .{ .geometry = .{ .sphere = .{ .radius = 0.25 } }, .material = .{ .restitution = 0.8 } });
    var highest_after: f32 = 0;
    var bounced = false;
    var t: f32 = 0;
    while (t < 2) : (t += dt) {
        try world.step(dt);
        const b = world.body(ball).?;
        if (b.linear_velocity.y > 0.5) bounced = true;
        if (bounced) highest_after = @max(highest_after, b.center.y);
    }
    try testing.expect(bounced);
    // Back up to most of the way.
    try testing.expect(highest_after > 1.5);
}

test "hulls and cylinders fall and rest on their faces" {
    var world = try withFloor();
    defer world.deinit();
    var corners: [8]Vec3 = undefined;
    for (&corners, 0..) |*v, i| v.* = .init(if (i & 1 != 0) 0.6 else -0.6, if (i & 2 != 0) 0.2 else -0.2, if (i & 4 != 0) 0.4 else -0.4);
    var slab = try Hull.init(testing.allocator, &corners);
    defer slab.deinit(testing.allocator);
    const rock = try world.createBody(.{ .position = .init(-2, 2, 0), .rotation = .fromAxisAngle(.init(1, 0, 0.3), 0.6) });
    _ = try world.addShape(rock, .{ .geometry = .{ .hull = &slab } });
    const drum = try world.createBody(.{ .position = .init(2, 2, 0) });
    _ = try world.addShape(drum, .cylinder(0.5, 0.4));
    try stepFor(&world, 4);
    try testing.expectApproxEqAbs(@as(f32, 0.2), world.body(rock).?.center.y, 0.03);
    try testing.expectApproxEqAbs(@as(f32, 0.5), world.body(drum).?.center.y, 0.03);
}

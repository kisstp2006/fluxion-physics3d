// SPDX-License-Identifier: BSD-2-Clause

//! Whole worlds held together by joints, stepped: pendulums that keep their
//! length, limits that hold, a ragdoll's chain falling onto a floor.

const std = @import("std");
const testing = std.testing;

const World = @import("World.zig");
const joint = @import("joint.zig");
const geometry = @import("geometry.zig");
const Vec3 = geometry.Vec3;
const Quat = geometry.Quat;
const Transform = geometry.Transform;

const dt = 1.0 / 60.0;

fn stepFor(world: *World, seconds: f32) !void {
    var t: f32 = 0;
    while (t < seconds) : (t += dt) try world.step(dt);
}

/// A ball a metre under a pin that never moves, pushed sideways, with a
/// joint of `kind` and the limits given.
fn pendulum(world: *World, def: joint.Def) !struct { pin: World.BodyId, bob: World.BodyId, joint: World.JointId } {
    const pin = try world.createBody(.{ .type = .static, .position = .init(0, 5, 0) });
    const bob = try world.createBody(.{ .position = .init(1, 5, 0), .linear_velocity = .init(0, 0, 3) });
    _ = try world.addShape(bob, .sphere(0.2));
    var made = def;
    made.body_a = pin;
    made.body_b = bob;
    return .{ .pin = pin, .bob = bob, .joint = try world.createJoint(made) };
}

test "a ball on a pin keeps its distance from the pin as it swings" {
    var world: World = .init(testing.allocator, .{});
    defer world.deinit();
    // The joint at the pin, seen from the bob a metre to its +x.
    const held = try pendulum(&world, .{ .body_a = .none, .body_b = .none, .frame_b = .at(.init(-1, 0, 0)) });
    var lowest: f32 = 5;
    for (0..240) |_| {
        try world.step(dt);
        const b = world.body(held.bob).?;
        try testing.expectApproxEqAbs(@as(f32, 1), b.center.sub(.init(0, 5, 0)).len(), 0.02);
        lowest = @min(lowest, b.center.y);
    }
    // It swung down under the pin.
    try testing.expect(lowest < 4.2);
}

test "a ball joint's cone stops a swing at its edge" {
    var world: World = .init(testing.allocator, .{});
    defer world.deinit();
    // Frames whose x points from the pin to the bob: at rest it sticks out
    // sideways, and it may fall no further than a fifth of a right angle.
    const cone = 0.2 * std.math.pi / 2.0;
    const held = try pendulum(&world, .{ .body_a = .none, .body_b = .none, .frame_b = .at(.init(-1, 0, 0)), .swing = cone });
    try stepFor(&world, 3);
    const j = world.joint(held.joint).?;
    const lean = j.swingAngle(world.body(held.pin).?, world.body(held.bob).?);
    try testing.expect(lean < cone + 0.03);
    // It leans on the edge, held there against gravity.
    try testing.expect(lean > cone - 0.05);
}

test "a ball joint's twist stops at either end" {
    var world: World = .init(testing.allocator, .{});
    defer world.deinit();
    const pin = try world.createBody(.{ .type = .static });
    const wheel = try world.createBody(.{ .position = .init(1, 0, 0), .gravity_scale = 0 });
    _ = try world.addShape(wheel, .box(.init(0.3, 0.3, 0.3)));
    const made = try world.createJoint(.{ .body_a = pin, .body_b = wheel, .frame_b = .at(.init(-1, 0, 0)), .twist = .{ .lower = -0.5, .upper = 0.25 } });
    for ([_]f32{ 6, -6 }) |spin| {
        world.body(wheel).?.setAngularVelocity(.init(spin, 0, 0));
        try stepFor(&world, 1);
        const turned = world.joint(made).?.twistAngle(world.body(pin).?, world.body(wheel).?);
        if (spin > 0) {
            try testing.expectApproxEqAbs(@as(f32, 0.25), turned, 0.04);
        } else {
            try testing.expectApproxEqAbs(@as(f32, -0.5), turned, 0.04);
        }
    }
}

test "a hinge turns only about its axis, and only as far as its limits" {
    var world: World = .init(testing.allocator, .{});
    defer world.deinit();
    // A door: the hinge's axis is z in its frame, turned to stand up along y.
    const upright = Quat.fromAxisAngle(.unit_x, -std.math.pi / 2.0);
    const frame = Transform{ .rotation = upright };
    const post = try world.createBody(.{ .type = .static });
    const door = try world.createBody(.{ .position = .init(0.5, 0, 0), .gravity_scale = 0 });
    _ = try world.addShape(door, .box(.init(0.5, 1, 0.05)));
    const made = try world.createJoint(.{
        .body_a = post,
        .body_b = door,
        .frame_a = frame,
        .frame_b = .{ .position = .init(-0.5, 0, 0), .rotation = upright },
        .kind = .hinge,
        .twist = .{ .lower = -1, .upper = 1 },
    });
    // Pushed every way at once: it swings about y, and no other way.
    world.body(door).?.setAngularVelocity(.init(2, 5, 2));
    try stepFor(&world, 2);
    const b = world.body(door).?;
    const up = b.transform.rotation.rotate(.unit_y);
    try testing.expect(up.dot(.unit_y) > 0.995);
    const turned = world.joint(made).?.twistAngle(world.body(post).?, b);
    try testing.expect(@abs(turned) < 1.04);
    // Its edge stays on the hinge.
    try testing.expect(b.transform.apply(.init(-0.5, 0, 0)).len() < 0.02);
}

test "two bodies a joint holds do not collide, unless it says so; a body taken away takes its joints" {
    var world: World = .init(testing.allocator, .{ .gravity = .zero });
    defer world.deinit();
    const a = try world.createBody(.{});
    const b = try world.createBody(.{ .position = .init(0.3, 0, 0) });
    _ = try world.addShape(a, .sphere(0.5));
    _ = try world.addShape(b, .sphere(0.5));
    const made = try world.createJoint(.{ .body_a = a, .body_b = b, .frame_b = .at(.init(-0.3, 0, 0)) });
    try stepFor(&world, 0.5);
    // Overlapping, and left so: no contact pushed them apart.
    try testing.expectEqual(@as(usize, 0), world.touchingCount());
    try testing.expectApproxEqAbs(@as(f32, 0.3), world.body(b).?.center.sub(world.body(a).?.center).len(), 0.01);
    try testing.expectError(error.SameBody, world.createJoint(.{ .body_a = a, .body_b = a }));
    world.destroyBody(b);
    try testing.expect(world.joint(made) == null);
    try testing.expectEqual(@as(usize, 0), world.jointCount());
    try testing.expectEqual(@as(usize, 0), world.joined.count());
}

test "a limp chain of capsules falls onto the floor, stays in one piece, and sleeps" {
    var world: World = .init(testing.allocator, .{});
    defer world.deinit();
    const ground = try world.createBody(.{ .type = .static, .position = .init(0, -0.5, 0) });
    _ = try world.addShape(ground, .box(.init(20, 0.5, 20)));
    // Five links lying along x, two metres up, each joined to the last at
    // its end, swinging up to half a right angle and twisting a little.
    const lying = Quat.fromAxisAngle(.unit_z, std.math.pi / 2.0);
    var links: [5]World.BodyId = undefined;
    for (&links, 0..) |*link, i| {
        link.* = try world.createBody(.{ .position = .init(@as(f32, @floatFromInt(i)) * 0.5, 2, 0), .angular_damping = 2 });
        _ = try world.addShape(link.*, .{ .geometry = .{ .capsule = .{ .half_height = 0.15, .radius = 0.08 } }, .offset = .{ .rotation = lying } });
        if (i > 0) _ = try world.createJoint(.{
            .body_a = links[i - 1],
            .body_b = link.*,
            .frame_a = .at(.init(0.25, 0, 0)),
            .frame_b = .at(.init(-0.25, 0, 0)),
            .swing = std.math.pi / 4.0,
            .twist = .{ .lower = -0.2, .upper = 0.2 },
        });
    }
    // Kicked, so it lands twisted.
    world.body(links[4]).?.setLinearVelocity(.init(0, 2, 3));
    try stepFor(&world, 6);
    for (links[1..], 1..) |link, i| {
        const before = world.body(links[i - 1]).?;
        const after = world.body(link).?;
        const gap = after.transform.apply(.init(-0.25, 0, 0)).sub(before.transform.apply(.init(0.25, 0, 0))).len();
        try testing.expect(gap < 0.02);
        try testing.expect(after.center.y > 0.05 and after.center.y < 0.2);
    }
    try testing.expectEqual(@as(usize, 0), world.awakeCount());
}

test "a joint started past its limit turns back into it no faster than the most" {
    var world: World = .init(testing.allocator, .{ .gravity = .zero });
    defer world.deinit();
    const pin = try world.createBody(.{ .type = .static });
    const arm = try world.createBody(.{ .position = .init(1, 0, 0), .rotation = .fromAxisAngle(.unit_x, 1.5) });
    _ = try world.addShape(arm, .box(.init(0.3, 0.1, 0.1)));
    const made = try world.createJoint(.{ .body_a = pin, .body_b = arm, .frame_b = .at(.init(-1, 0, 0)), .twist = .{ .lower = -0.2, .upper = 0.2 } });
    try world.step(dt);
    try testing.expect(world.body(arm).?.angular_velocity.len() <= world.settings.max_limit_turn + 0.01);
    try stepFor(&world, 1);
    const turned = world.joint(made).?.twistAngle(world.body(pin).?, world.body(arm).?);
    try testing.expectApproxEqAbs(@as(f32, 0.2), turned, 0.05);
}

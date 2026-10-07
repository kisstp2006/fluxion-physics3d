// SPDX-License-Identifier: BSD-2-Clause

//! Fluxion Physics 3D - rigid bodies in space.
//!
//!   `World`     the bodies, the shapes on them, and one step of time
//!   `Body`      a rigid body: where it is, how it moves, what it weighs
//!   `shape`     balls, boxes, capsules, cylinders, hulls and meshes
//!   `Hull`      a convex solid built from any points
//!   `Mesh`      triangles to collide with, for what never moves
//!   `collide`   where two solids touch, and how deep
//!   `gjk`       how near two convex solids come, and where a cast stops
//!   `contact`   a contact as the solver sees it
//!   `Tree`      a tree of boxes, for finding what is near
//!   `geometry`  places, turns, boxes, and the matrix sums
//!
//! ```zig
//! const physics3d = @import("fluxion_physics3d");
//!
//! var world: physics3d.World = .init(gpa, .{});
//! defer world.deinit();
//! const ground = try world.createBody(.{ .type = .static });
//! _ = try world.addShape(ground, .box(.init(10, 0.5, 10)));
//! const ball = try world.createBody(.{ .position = .init(0, 5, 0) });
//! _ = try world.addShape(ball, .sphere(0.5));
//! try world.step(1.0 / 60.0);
//! const where = world.body(ball).?.position();
//! ```
//!
//! **It knows nothing about entities.** A body is named by a handle, and a
//! game keeps the handle wherever it keeps things. That is what lets it be
//! built and tested on its own, and used by something that is not the
//! engine.

const std = @import("std");

pub const World = @import("World.zig");
pub const Body = @import("body.zig");
pub const shape = @import("shape.zig");
pub const Hull = @import("hull.zig");
pub const Mesh = @import("mesh.zig");
pub const collide = @import("collide.zig");
pub const gjk = @import("gjk.zig");
pub const contact = @import("contact.zig");
pub const Tree = @import("tree.zig");
pub const geometry = @import("geometry.zig");
pub const Softness = @import("Softness.zig");

pub const BodyId = Body.Id;
pub const ShapeId = Body.ShapeId;
pub const BodyType = Body.Type;
pub const BodyDef = Body.Def;

pub const Shape = shape.Shape;
pub const Geometry = shape.Geometry;
pub const Material = shape.Material;
pub const Filter = shape.Filter;
pub const FilterRule = shape.FilterRule;
pub const Mix = shape.Mix;

pub const Settings = World.Settings;
pub const ContactEvent = World.ContactEvent;
pub const RayHit = World.RayHit;
pub const ShapeHit = World.ShapeHit;
pub const QueryFilter = World.QueryFilter;
pub const Penetration = World.Penetration;
pub const Manifold = collide.Manifold;

pub const Vec3 = geometry.Vec3;
pub const Quat = geometry.Quat;
pub const Transform = geometry.Transform;
pub const Aabb = geometry.Aabb;

/// The arithmetic underneath.
pub const math = @import("fluxion_math");

test {
    _ = World;
    _ = Body;
    _ = shape;
    _ = Hull;
    _ = Mesh;
    _ = collide;
    _ = gjk;
    _ = contact;
    _ = Tree;
    _ = geometry;
    _ = Softness;
    _ = @import("world_test.zig");
}

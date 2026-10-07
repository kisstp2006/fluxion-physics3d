# Fluxion Physics 3D

Rigid bodies in space. For Zig 0.16, on any machine and in a browser.

| Module | What it is |
| --- | --- |
| `World` | The bodies, the shapes on them, one step of time, and what to ask of it. |
| `Body` | A rigid body: where it is, how it moves, what it weighs, whether it sleeps. |
| `shape` | Balls, boxes, capsules, cylinders, hulls and meshes; materials and filters. |
| `Hull` | A convex solid built from any points. |
| `Mesh` | Triangles to collide with, for what never moves. |
| `collide` | Where two solids touch, and how deep. |
| `gjk` | How near two convex solids come, and where a cast of one stops. |
| `contact` | A contact as the solver sees it. |
| `Softness` | A spring as the solver steps it, and why every contact is a little of one. |
| `Tree` | A tree of boxes: what is near, without looking at all of it. |
| `geometry` | Places, turns, boxes, and the few matrix sums the rest needs. |

```zig
const physics3d = @import("fluxion_physics3d");

var world: physics3d.World = .init(gpa, .{});
defer world.deinit();

const ground = try world.createBody(.{ .type = .static, .position = .init(0, -0.5, 0) });
_ = try world.addShape(ground, .box(.init(20, 0.5, 20)));

const crate = try world.createBody(.{ .position = .init(0, 4, 0) });
_ = try world.addShape(crate, .{ .geometry = .{ .box = .{ .half = .init(0.5, 0.5, 0.5) } }, .material = .{ .friction = 0.8 } });

const player = try world.createBody(.{ .position = .init(2, 1, 0), .lock_rotation = true });
_ = try world.addShape(player, .capsule(0.6, 0.3));

try world.step(1.0 / 60.0);
const where = world.body(crate).?.position();
const below = world.castRay(where, .init(0, -10, 0), .{ .ignore = crate });
```

**A body is a handle, and so is a shape.** Both are
[Fluxion Id](https://github.com/kisstp2006/fluxion-id) generational handles:
eight bytes to copy into a component, and a handle to something destroyed
answers null for ever, however many times its slot is reused.
`world.body(h)` hands out a pointer to push on; hold the handle, not the
pointer.

**Three kinds of body.** Static never moves and holds everything up.
Kinematic moves the way it is told and is pushed by nothing: a lift, a door.
Dynamic is what physics is for. Two bodies that cannot both move never
collide - a sensor sees them all - and a kinematic platform carries what
stands on it.

**Six geometries.** A ball, a box, a capsule and a cylinder - the last two
standing along their own `y`, as a character stands - a convex hull built
from any points (`Hull.init` keeps the ones on the outside and makes the
triangles that lie in one plane one face, so a hull of a box's corners has
six faces), and a mesh of triangles for what never moves: the level. A
cylinder collides as a prism of sixteen sides. Several shapes on one body,
each placed by its `offset`, make one solid.

**Mass comes from density and volume**, so a big crate weighs more than a
small one without anyone typing a mass, and a body's centre of mass and
inertia are where and what its shapes say. A body moves about its centre of
mass; its origin - where the game put it - follows.

**The axes are the engine's.** `+y` up, right-handed, metres, kilograms
and seconds; gravity pulls towards `-y`.

## Touching

**Two solids touch at up to four points**, each with how far apart the
surfaces are there, along one normal. Two round solids - balls and capsules,
each a point or a segment and a radius - meet at the nearest points of
their cores, and two capsules lying side by side at both ends of where they
overlap. A round solid against one with corners is the nearest points of its
core and the solid, by GJK, and a capsule lying along a face touches it at
both ends of the part over the face. Two solids with corners are tried along
every way they could be apart - each face of either, each edge of one
crossed with each edge of the other - and pushed apart along the way they
are least sunk; along a face, the face of the other that faces it most is
clipped to it, and the four of its corners that span the most are kept.

**Points are found a little before the surfaces meet**
(`speculative_distance`), so a falling body is stopped at the floor rather
than found already through it.

**A mesh is each of its triangles**, found through a tree of its own, and
solid from both sides. Only a static body may have one: a mesh has no
inside, so no mass.

## The step

**Soft and in substeps.** A step is cut into a few substeps, each
integrating the velocities, solving every contact once, moving the bodies
and relaxing once; how deep each point has sunk is worked out again every
substep from where the bodies now are. Sinking in is taken back as a very
stiff, damped spring would take it - Erin Catto's soft step; see `Softness`
- and the relaxing pass takes the push back out of the velocities, so a box
made half inside another slides out and stops. Friction pulls both ways
along the surface, no harder in all than the friction times the push.
Restitution comes last, from the speed each point came in at.

**Bodies fall asleep** in islands: everything touching everything else
that touches it, still for `time_to_sleep`. A step passes a sleeper over
until something wakes it - a push, a new speed, a move, or a body running
into it. A body riding something that moves never sleeps.

**What began and stopped touching** is said after each step:
`beginEvents` and `endEvents`, sensors included, a mesh's triangles counted
once.

## Asking

| | |
| --- | --- |
| `castRay(origin, translation, filter)` | The first shape a ray meets, where, and the normal there. |
| `castShape(geometry, transform, translation, filter)` | Where a shape moved along a line first comes to touch something, stopping the slop short. |
| `penetrations(geometry, transform, filter, out)` | How deep a shape has sunk into each it overlaps, and the way out: what a character steps out of walls with. |
| `overlapShape`, `overlapPoint`, `overlapAabb` | What a shape, a point or a box overlaps. |

A `QueryFilter` says which layers are seen, a body to pass over - the one
cast from - and whether sensors count.

## Tests

`zig build test` steps whole worlds - a box dropped and asleep, a stack of
six, a ball rolling down a slope, a capsule standing on a floor of
triangles, a ball through a sensor, a box riding a lift, a bouncing ball -
and builds the library for `wasm32-freestanding` too.

## Licence

BSD-2-Clause. See `LICENSE`.

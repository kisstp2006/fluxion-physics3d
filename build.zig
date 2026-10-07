// SPDX-License-Identifier: BSD-2-Clause

const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = module(b, target, optimize, .exported);

    // zig build test
    const tests = b.addTest(.{ .name = "fluxion-physics3d-tests", .root_module = mod });
    const test_step = b.step("test", "Run the library test suite");
    test_step.dependOn(&b.addRunArtifact(tests).step);

    // The same library for a browser, built as part of the tests: a change
    // that breaks the no-thread, no-operating-system path should fail here,
    // not in a game's web export later.
    const wasm_target = b.resolveTargetQuery(.{ .cpu_arch = .wasm32, .os_tag = .freestanding });
    const wasm = b.addLibrary(.{
        .name = "fluxion-physics3d-wasm",
        .root_module = module(b, wasm_target, .ReleaseSmall, .private),
    });
    test_step.dependOn(&wasm.step);

    // zig build docs -> zig-out/docs
    const docs_lib = b.addLibrary(.{ .name = "fluxion-physics3d", .root_module = mod });
    b.step("docs", "Generate API documentation into zig-out/docs").dependOn(&b.addInstallDirectory(.{
        .source_dir = docs_lib.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "docs",
    }).step);
}

/// The importable module, for one target. Consumers do:
///   const physics3d = @import("fluxion_physics3d");
fn module(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    /// `.private` for the browser build, which is made for wasm whatever the
    /// target and must not take the name a consumer asks for.
    kind: enum { exported, private },
) *std.Build.Module {
    const math = b.dependency("fluxion_math", .{ .target = target, .optimize = optimize });
    const ident = b.dependency("fluxion_id", .{ .target = target, .optimize = optimize });

    const options: std.Build.Module.CreateOptions = .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "fluxion_math", .module = math.module("fluxion_math") },
            .{ .name = "fluxion_id", .module = ident.module("fluxion_id") },
        },
    };
    return if (kind == .exported) b.addModule("fluxion_physics3d", options) else b.createModule(options);
}

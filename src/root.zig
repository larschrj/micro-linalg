const std = @import("std");

pub const matrix = @import("matrix.zig");
pub const lu = @import("lu.zig");

pub const DenseMatrix = matrix.DenseMatrix;
pub const LuSolver = lu.LuSolver;

test {
    std.testing.refAllDecls(@This());
}

const std = @import("std");

pub fn DenseMatrix(comptime N: usize, comptime M: usize, comptime T: type) type {
    return struct {
        data: [N][M]T,

        const Self = @This();

        pub fn init(data: [N][M]T) Self {
            return .{ .data = data };
        }

        pub fn initZeros() Self {
            return .{ .data = std.mem.zeroes([N][M]T) };
        }

        pub fn scalarMul(self: *Self, scalar: T) void {
            for (0..M) |j| {
                for (0..N) |i| {
                    self.data[i][j] *= scalar;
                }
            }
        }
    };
}

const std = @import("std");

pub const LuError = error{
    SingularMatrix,
};

pub fn LuSolver(comptime N: usize, comptime T: type) type {
    return struct {
        lu: [N][N]T,
        p: [N]usize,

        const Self = @This();

        pub fn init(input: [N][N]T) Self {
            return .{
                .lu = input,
                .p = undefined,
            };
        }

        // LU decomposition with row pivoting
        pub fn decompose(self: *Self) LuError!void {
            // Pivot rows are initialised to the existing row order
            for (0..N) |i| {
                self.p[i] = i;
            }

            for (0..N) |k| {
                var max_value: T = 0.0;
                var pivot_index: usize = self.p[k];

                // 1. Find pivot row
                for (k..N) |i| {
                    const value = @abs(self.lu[self.p[i]][k]);

                    if (value > max_value) {
                        pivot_index = i;
                        max_value = value;
                    }
                }

                if (max_value == 0.0) {
                    return LuError.SingularMatrix;
                }

                // 2. Swap permutation indices
                const temp = self.p[k];
                self.p[k] = self.p[pivot_index];
                self.p[pivot_index] = temp;

                // 3. Perform decomposition step
                const pivot_row = self.p[k];
                for ((k + 1)..N) |i| {
                    const current_row = self.p[i];

                    // Row scale factor
                    self.lu[current_row][k] /= self.lu[pivot_row][k];

                    // Gaussian elimination (L and U)
                    for ((k + 1)..N) |j| {
                        self.lu[current_row][j] -= self.lu[current_row][k] * self.lu[pivot_row][j];
                    }
                }
            }
        }

        pub fn solve(self: *Self, b: [N]T) [N]T {
            var x: [N]T = .{0} ** N;

            // L*y = b
            for (0..N) |i| {
                const current_row = self.p[i];
                x[i] = b[current_row];

                for (0..i) |j| {
                    x[i] -= self.lu[current_row][j] * x[j];
                }
            }

            // U*x = y
            for (0..N) |i| {
                const current_row = self.p[N - i - 1];
                x[N - i - 1] = x[N - i - 1];

                for (0..i) |j| {
                    x[N - i - 1] -= self.lu[current_row][N - j - 1] * x[N - j - 1];
                }

                x[N - i - 1] /= self.lu[current_row][N - i - 1];
            }

            return x;
        }
    };
}

test "LU test" {
    const A = [3][3]f16{ [_]f16{ 2, 7, 5 }, [_]f16{ -2, 3, -2 }, [_]f16{ 4, -2, 2 } };

    const LU3_f16 = LuSolver(3, f16);
    var lu3_f16 = LU3_f16.init(A);
    try lu3_f16.decompose();

    std.debug.print("\nLU: {any}", .{lu3_f16.lu});
    std.debug.print("\nP: {any}", .{lu3_f16.p});

    var b = [_]f16{ 14, -1, 4 };
    std.debug.print("\nb: {any}", .{b});
    var x = lu3_f16.solve(b);
    std.debug.print("\nx: {any}", .{x});

    b = [_]f16{ 12, -13, 16 };
    std.debug.print("\nb: {any}", .{b});
    x = lu3_f16.solve(b);
    std.debug.print("\nx: {any}", .{x});

    const B = [3][3]f32{ [_]f32{ 1, 2, 3 }, [_]f32{ 4, 5, 6 }, [_]f32{ 7, 8, 9.01 } };
    const LU3_f32 = LuSolver(3, f32);
    var lu3_f32 = LU3_f32.init(B);
    try lu3_f32.decompose();
    std.debug.print("\nLU: {any}", .{lu3_f32.lu});
    std.debug.print("\nP: {any}", .{lu3_f32.p});

    const c = [_]f32{ 6, 15, 24.01 };
    std.debug.print("\nb: {any}", .{c});
    const y = lu3_f32.solve(c);
    std.debug.print("\nx: {any}", .{y});
}

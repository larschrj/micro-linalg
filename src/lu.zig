const std = @import("std");

pub const LuError = error{
    SingularMatrix,
    IllConditioned,
};

pub fn LuSolver(comptime N: usize, comptime T: type) type {
    return struct {
        lu: [N][N]T,
        p: [N]usize,
        pivot_ratio: T,
        max_pivot_ratio: T = 1.0 / @sqrt(std.math.floatEps(T)),

        const Self = @This();

        pub fn init(input: [N][N]T) Self {
            return .{
                .lu = input,
                .p = undefined,
                .pivot_ratio = 1.0,
            };
        }

        // LU decomposition with row pivoting
        pub fn decompose(self: *Self) LuError!void {
            var max_pivot: T = 0.0;
            var min_pivot: T = 0.0;

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

                if (k == 0) {
                    max_pivot = max_value;
                    min_pivot = max_value;
                } else {
                    if (max_value > max_pivot) max_pivot = max_value;
                    if (max_value < min_pivot) min_pivot = max_value;
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

            self.pivot_ratio = max_pivot / min_pivot;

            if (self.pivot_ratio > self.max_pivot_ratio) {
                return LuError.IllConditioned;
            }
        }

        pub fn solve(self: *Self, b: *const [N]T, x: *[N]T) void {
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
        }

        pub fn solveMatrix(self: *Self, comptime M: usize, B: *const [N][M]T, X: *[N][M]T) void {
            // Check B has N rows
            if (B.*.len != N) {
                @compileError("Matrix B row count must match LU row/column count");
            }

            for (0..M) |j| {
                var b_col: [N]T = undefined;
                for (0..N) |i| {
                    b_col[i] = B[i][j];
                }

                var x_col: [N]T = undefined;
                self.solve(&b_col, &x_col);
                for (0..N) |i| {
                    X[i][j] = x_col[i];
                }
            }
        }
    };
}

fn expectApproxEqSlice(comptime T: type, expected: []const T, actual: []const T, tolerance: T) !void {
    try std.testing.expectEqual(expected.len, actual.len);
    for (expected, 0..) |e, i| {
        try std.testing.expectApproxEqAbs(e, actual[i], tolerance);
    }
}

test "LuSolver - decompose and solve well-conditioned f16 matrix" {
    const A = [3][3]f16{ [_]f16{ 2, 7, 5 }, [_]f16{ -2, 3, -2 }, [_]f16{ 4, -2, 2 } };
    var solver = LuSolver(3, f16).init(A);
    try solver.decompose();

    const b1 = [_]f16{ 14, -1, 4 };
    var x1: [3]f16 = undefined;
    solver.solve(&b1, &x1);
    try expectApproxEqSlice(f16, &[_]f16{ 1, 1, 1 }, &x1, 0.001);

    const b2 = [_]f16{ 12, -13, 16 };
    var x2: [3]f16 = undefined;
    solver.solve(&b2, &x2);
    try expectApproxEqSlice(f16, &[_]f16{ 2, -1, 3 }, &x2, 0.001);
}

test "LuSolver - decompose and solve well-conditioned f32 matrix" {
    const B = [3][3]f32{ [_]f32{ 1, 2, 3 }, [_]f32{ 4, 5, 6 }, [_]f32{ 7, 8, 9.01 } };
    var solver = LuSolver(3, f32).init(B);
    try solver.decompose();

    const c = [_]f32{ 6, 15, 24.01 };
    var y: [3]f32 = undefined;
    solver.solve(&c, &y);

    try expectApproxEqSlice(f32, &[_]f32{ 1, 1, 1 }, &y, 0.001);
}

test "LuSolver - ill-conditioned f32 matrix returns error" {
    const B = [3][3]f32{ [_]f32{ 1, 2, 3 }, [_]f32{ 4, 5, 6 }, [_]f32{ 7, 8, 9.001 } };
    var solver = LuSolver(3, f32).init(B);

    const err = solver.decompose();
    try std.testing.expectError(LuError.IllConditioned, err);
}

test "LuSolver - solve matrix" {
    const A = [3][3]f32{ [_]f32{ 1, 2, 3 }, [_]f32{ 4, 5, 6 }, [_]f32{ 7, 8, 9.01 } };
    var solver = LuSolver(3, f32).init(A);
    try solver.decompose();

    const B = [3][2]f32{
        [_]f32{ 6, 12 },
        [_]f32{ 15, 30 },
        [_]f32{ 24.01, 48.02 },
    };
    var X: [3][2]f32 = undefined;

    solver.solveMatrix(2, &B, &X);

    try expectApproxEqSlice(f32, &[_]f32{ 1, 2 }, &X[0], 0.001);
    try expectApproxEqSlice(f32, &[_]f32{ 1, 2 }, &X[1], 0.001);
    try expectApproxEqSlice(f32, &[_]f32{ 1, 2 }, &X[2], 0.001);
}

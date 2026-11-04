const std = @import("std");
const SymmetricMatrix = @import("matrix.zig").SymmetricMatrix;

pub const CholeskyError = error{
    NotPositiveDefinite,
};

pub fn CholeskySolver(comptime N: usize, comptime T: type) type {
    return struct {
        L: SymmetricMatrix(N, T),

        const Self = @This();

        pub fn init(A: SymmetricMatrix(N, T)) Self {
            return .{ .L = A };
        }

        pub fn decompose(self: *Self) CholeskyError!void {
            for (0..N) |i| {
                var diag_sum: T = 0;

                for (0..i) |j| {
                    var sum: T = 0.0;
                    for (0..j) |k| {
                        sum += self.L.get(i, k) * self.L.get(j, k);
                    }
                    const Lij = (self.L.get(i, j) - sum) / self.L.get(j, j);
                    self.L.set(i, j, Lij);
                    diag_sum += Lij * Lij;
                }

                const Lii_sq: T = self.L.get(i, i) - diag_sum;
                if (Lii_sq <= 0) {
                    return CholeskyError.NotPositiveDefinite;
                }
                self.L.set(i, i, @sqrt(Lii_sq));
            }
        }
    };
}

test "Cholesky" {
    const Sym3f32 = SymmetricMatrix(3, f32);
    const Solver3f32 = CholeskySolver(3, f32);

    const A: Sym3f32 = .{ .data = [_]f32{ 9.0, 3.0, 10.0, 3.0, 4.0, 11.0 } };
    var solver = Solver3f32.init(A);
    try solver.decompose();

    const expected: Sym3f32 = .{ .data = [_]f32{ 3.0, 1.0, 3.0, 1.0, 1.0, 3.0 } };
    for (expected.data, solver.L.data) |exp, actual| {
        try std.testing.expectApproxEqAbs(exp, actual, std.math.floatEps(f32));
    }
}

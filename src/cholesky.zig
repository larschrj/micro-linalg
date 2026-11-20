const std = @import("std");
const SymmetricMatrix = @import("matrix.zig").SymmetricMatrix;

pub const CholeskyError = error{
    NotPositiveDefinite,
};

pub fn CholeskySolver(comptime N: usize, comptime T: type) type {
    return struct {
        L: SymmetricMatrix(N, T),

        const Self = @This();

        pub fn decompose(A: SymmetricMatrix(N, T)) CholeskyError!Self {
            var L = A;

            for (0..N) |i| {
                var diag_sum: T = 0;

                for (0..i) |j| {
                    var sum: T = 0.0;
                    for (0..j) |k| {
                        sum += L.get(i, k) * L.get(j, k);
                    }
                    const Lij = (L.get(i, j) - sum) / L.get(j, j);
                    L.set(i, j, Lij);
                    diag_sum += Lij * Lij;
                }

                const Lii_sq: T = L.get(i, i) - diag_sum;
                // Catch negative and NaN values
                if (!(Lii_sq > 0)) {
                    return CholeskyError.NotPositiveDefinite;
                }
                L.set(i, i, @sqrt(Lii_sq));
            }

            return .{ .L = L };
        }
    };
}

test "Cholesky decomposition of symmetric positive definite matrix" {
    const Sym3f32 = SymmetricMatrix(3, f32);
    const Solver3f32 = CholeskySolver(3, f32);

    const A: Sym3f32 = .{ .data = [_]f32{ 9.0, 3.0, 10.0, 3.0, 4.0, 11.0 } };
    const solver = try Solver3f32.decompose(A);

    const expected: Sym3f32 = .{ .data = [_]f32{ 3.0, 1.0, 3.0, 1.0, 1.0, 3.0 } };
    for (expected.data, solver.L.data) |exp, actual| {
        try std.testing.expectApproxEqAbs(exp, actual, std.math.floatEps(f32));
    }
}

test "Cholesky decomposition fails for symmetric non-positive definite matrix" {
    const Sym2f32 = SymmetricMatrix(2, f32);
    const Solver2f32 = CholeskySolver(2, f32);

    // Indefinite matrix: [1, 2; 2, 1]
    const A: Sym2f32 = .{ .data = [_]f32{ 1.0, 2.0, 1.0 } };

    try std.testing.expectError(CholeskyError.NotPositiveDefinite, Solver2f32.decompose(A));
}

test "Cholesky decompostion fails for symmetric semi-definite matrix" {
    const Sym2f32 = SymmetricMatrix(2, f32);
    const Solver2f32 = CholeskySolver(2, f32);

    // Indefinite matrix: [1, 1; 1, 1]
    const A: Sym2f32 = .{ .data = [_]f32{ 1.0, 1.0, 1.0 } };

    try std.testing.expectError(CholeskyError.NotPositiveDefinite, Solver2f32.decompose(A));
}

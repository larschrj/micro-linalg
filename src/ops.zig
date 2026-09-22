const std = @import("std");
const matrix = @import("matrix.zig");
const DenseMatrix = matrix.DenseMatrix;
const SymmetricMatrix = matrix.SymmetricMatrix;

/// Multiplies the matrices `lhs` and `rhs` and stores the result in `out`.
/// The number of columns in `lhs` must equal the number of rows in `rhs`.
/// `lhs` and `rhs` must be a pointer to a matrix.
/// `out` must be a pointer to a DenseMatrix.
pub fn multiplyMatrixMatrix(lhs: anytype, rhs: anytype, out: anytype) void {
    const LhsType = @TypeOf(lhs);
    const RhsType = @TypeOf(rhs);
    const OutType = @TypeOf(out);

    // --- Type checks ---
    // lhs
    if (@typeInfo(LhsType) != .pointer) {
        @compileError("Expected 'lhs' to be a pointer to matrix, but received a value of type'" ++ @typeName(LhsType) ++ "'. Did you forget to pass it by reference using '&'?");
    }
    const LhsChildType = @typeInfo(LhsType).pointer.child;
    const lhs_is_dense_or_sym = @hasDecl(LhsChildType, "is_dense_matrix") or @hasDecl(LhsChildType, "is_symmetric_matrix");
    if (!lhs_is_dense_or_sym) {
        @compileError("Expected 'lhs' to be a pointer to DenseMatrix or SymmetricMatrix, but received a pointer to type'" ++ @typeName(LhsType) ++ "'.");
    }

    // rhs
    if (@typeInfo(RhsType) != .pointer) {
        @compileError("Expected 'rhs' to be a pointer to matrix, but received a value of type'" ++ @typeName(RhsType) ++ "'. Did you forget to pass it by reference using '&'?");
    }
    const RhsChildType = @typeInfo(RhsType).pointer.child;
    const rhs_is_dense_or_sym = @hasDecl(RhsChildType, "is_dense_matrix") or @hasDecl(RhsChildType, "is_symmetric_matrix");
    if (!rhs_is_dense_or_sym) {
        @compileError("Expected 'rhs' to be a pointer to DenseMatrix or SymmetricMatrix, but received a pointer to type'" ++ @typeName(RhsType) ++ "'.");
    }

    // out
    if (@typeInfo(OutType) != .pointer) {
        @compileError("Expected 'out' to be a pointer to DenseMatrix, but received a value of type'" ++ @typeName(OutType) ++ "'. Did you forget to pass it by reference using '&'?");
    }
    const OutChildType = @typeInfo(OutType).pointer.child;
    if (!@hasDecl(OutChildType, "is_dense_matrix")) {
        @compileError("Expected 'out' to be a pointer to DenseMatrix, but received a pointer to type'" ++ @typeName(OutChildType) ++ "'.");
    }

    // --- Dimension checks ---
    const ln = @TypeOf(lhs.*).rows;
    const lm = @TypeOf(lhs.*).columns; // lhs columns must match rhs rows
    const rn = @TypeOf(rhs.*).rows;
    const rm = @TypeOf(rhs.*).columns;
    const on = @TypeOf(out.*).rows;
    const om = @TypeOf(out.*).columns;
    if (lm != rn) {
        @compileError("lhs column dimension not equal to rhs row dimension");
    }
    if (on != ln) {
        @compileError("out row dimension not equal to lhs row dimension");
    }
    if (om != rm) {
        @compileError("out column dimension not equal to rhs column dimension");
    }

    // --- Calculate out row by row, using lhs to scale and add rows of rhs ---
    for (0..ln) |i| {
        // 1. Initialize row i
        const a_i0 = lhs.get(i, 0);
        for (0..rm) |j| {
            out.set(i, j, a_i0 * rhs.get(0, j));
        }

        // 2. Accumulate result in row i
        for (1..lm) |k| {
            const a_ik = lhs.get(i, k);
            for (0..rm) |j| {
                const current = out.*.get(i, j);
                out.set(i, j, current + a_ik * rhs.get(k, j));
            }
        }
    }
}

test "Multiply DenseMatrix with DenseMatrix" {
    const d3x3f32 = DenseMatrix(3, 3, f32);
    const mat_a: d3x3f32 = .{ .data = [3][3]f32{
        [_]f32{ 1.0, 2.0, 3.0 },
        [_]f32{ 4.0, 5.0, 6.0 },
        [_]f32{ 7.0, 8.0, 9.0 },
    } };
    const mat_b: d3x3f32 = .{ .data = [3][3]f32{
        [_]f32{ 2.0, 3.0, 4.0 },
        [_]f32{ 5.0, 6.0, 7.0 },
        [_]f32{ 8.0, 9.0, 10.0 },
    } };
    var mat_c: d3x3f32 = .{ .data = [3][3]f32{
        [_]f32{ 0.0, 0.0, 0.0 },
        [_]f32{ 0.0, 0.0, 0.0 },
        [_]f32{ 0.0, 0.0, 0.0 },
    } };
    var mat_c_expected: d3x3f32 = .{ .data = [3][3]f32{
        [_]f32{ 36.0, 42.0, 48.0 },
        [_]f32{ 81.0, 96.0, 111.0 },
        [_]f32{ 126.0, 150.0, 174.0 },
    } };
    multiplyMatrixMatrix(&mat_a, &mat_b, &mat_c);

    for (0..d3x3f32.rows) |i| {
        for (0..d3x3f32.columns) |j| {
            try std.testing.expectApproxEqRel(mat_c.get(i, j), mat_c_expected.get(i, j), std.math.floatEps(f32));
        }
    }
}

const std = @import("std");

pub fn DenseMatrix(comptime N: usize, comptime M: usize, comptime T: type) type {
    return struct {
        data: [N][M]T,

        pub const is_dense_matrix = true;
        pub const rows = N;
        pub const columns = M;
        pub const Scalar = T;

        const Self = @This();

        pub fn init(data: [N][M]T) Self {
            return .{ .data = data };
        }

        pub fn initZeros() Self {
            return .{ .data = std.mem.zeroes([N][M]T) };
        }

        pub fn get(self: *const Self, row_index: usize, column_index: usize) T {
            return self.data[row_index][column_index];
        }

        pub fn set(self: *Self, row_index: usize, column_index: usize, value: T) void {
            self.data[row_index][column_index] = value;
        }
    };
}

pub fn SymmetricMatrix(comptime N: usize, comptime T: type) type {
    const packed_len = (N * (N + 1)) / 2;

    return struct {
        data: [packed_len]T,

        pub const is_symmetric_matrix = true;
        pub const rows = N;
        pub const columns = N;
        pub const Scalar = T;

        const Self = @This();

        pub fn init(data: [packed_len]T) Self {
            return .{ .data = data };
        }

        pub fn initZeros() Self {
            return .{ .data = std.mem.zeroes([packed_len]T) };
        }

        pub fn indexOf(row_index: usize, column_index: usize) usize {
            // Check bounds
            std.debug.assert(row_index < N);
            std.debug.assert(column_index < N);

            // Upper half of matrix is the same as lower half with indices swapped
            const row = @max(row_index, column_index);
            const col = @min(row_index, column_index);

            // Index of first element in row with column offset
            return (row * (row + 1)) / 2 + col;
        }

        pub fn get(self: *const Self, row_index: usize, column_index: usize) T {
            return self.data[indexOf(row_index, column_index)];
        }

        pub fn set(self: *Self, row_index: usize, column_index: usize, value: T) void {
            self.data[indexOf(row_index, column_index)] = value;
        }
    };
}

pub fn LowerTriangularMatrix(comptime N: usize, comptime T: type) type {
    const packed_len = (N * (N + 1)) / 2;

    return struct {
        data: [packed_len]T,

        pub const is_lower_triangular_matrix = true;
        pub const rows = N;
        pub const columns = N;
        pub const Scalar = T;

        const Self = @This();

        pub fn init(data: [packed_len]T) Self {
            return .{ .data = data };
        }

        pub fn initZeros() Self {
            return .{ .data = std.mem.zeroes([packed_len]T) };
        }

        pub fn indexOf(row_index: usize, column_index: usize) usize {
            std.debug.assert(row_index >= column_index);
            std.debug.assert(row_index < N);
            return (row_index * (row_index + 1)) / 2 + column_index;
        }

        pub fn get(self: *const Self, row_index: usize, column_index: usize) T {
            std.debug.assert(row_index < N);
            std.debug.assert(column_index < N);
            if (column_index > row_index) {
                return 0;
            } else {
                return self.data[indexOf(row_index, column_index)];
            }
        }

        pub fn set(self: *Self, row_index: usize, column_index: usize, value: T) void {
            self.data[indexOf(row_index, column_index)] = value;
        }
    };
}

test "Construct SymmetricMatrix" {
    const S3f32 = SymmetricMatrix(3, f32);
    var s0 = S3f32.init(.{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });

    // Get diagonal
    try std.testing.expectApproxEqAbs(5.0, s0.get(2, 2), std.math.floatEps(f32));

    // Set diagonal
    s0.set(2, 2, 6.0);
    try std.testing.expectApproxEqAbs(6.0, s0.get(2, 2), std.math.floatEps(f32));

    // Symmetric set
    s0.set(0, 2, 7.0);
    try std.testing.expectApproxEqAbs(7.0, s0.get(0, 2), std.math.floatEps(f32));
    try std.testing.expectApproxEqAbs(7.0, s0.get(2, 0), std.math.floatEps(f32));
}

test "SymmetricMatrix bug test" {
    const S3f32 = SymmetricMatrix(3, f32);
    var s0 = S3f32.init(.{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });

    // Setting and getting across the diagonal
    s0.set(0, 1, 99.0);
    try std.testing.expectApproxEqAbs(s0.get(1, 0), 99.0, std.math.floatEps(f32));
}

test "LowerTriangularMatrix" {
    const L3f32 = LowerTriangularMatrix(3, f32);
    var l0 = L3f32.init((.{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 }));

    // Get diagonal
    try std.testing.expectApproxEqAbs(5.0, l0.get(2, 2), std.math.floatEps(f32));

    // Get 0 above diagonal
    try std.testing.expectApproxEqAbs(0.0, l0.get(0, 2), std.math.floatEps(f32));

    // Set diagonal
    l0.set(0, 0, 6.0);
    try std.testing.expectApproxEqAbs(6.0, l0.get(0, 0), std.math.floatEps(f32));

    // Set below diagonal
    l0.set(2, 0, 7.0);
    try std.testing.expectApproxEqAbs(7.0, l0.get(2, 0), std.math.floatEps(f32));
}

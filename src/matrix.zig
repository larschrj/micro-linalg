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
            // Upper half of matrix is the same as lower half with indices swapped
            const swap_indices = column_index > row_index;
            const row = if (swap_indices) column_index else row_index;
            const col = if (swap_indices) row_index else column_index;

            // Row_index 1 is row number 2
            const row_count = row + 1;
            // Index of diagonal element in row with row_index
            const diagonal_index = @divFloor(row_count * (row_count + 1), 2) - 1;
            const index = if (row == col) {
                // Diagonal element
                return diagonal_index;
            } else {
                // Off diagonal element
                const offset = row - col;
                // Count backwards from diagonal element in row
                return diagonal_index - offset;
            };

            return index;
        }

        pub fn get(self: *const Self, row_index: usize, column_index: usize) T {
            return self.data[indexOf(row_index, column_index)];
        }

        pub fn set(self: *Self, row_index: usize, column_index: usize, value: T) void {
            self.data[indexOf(row_index, column_index)] = value;
        }
    };
}

test "Construct SymmetricMatrix" {
    const S3f32 = SymmetricMatrix(3, f32);
    var s0 = S3f32.init(.{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
    try std.testing.expectApproxEqAbs(s0.get(2, 2), 5.0, std.math.floatEps(f32));
    s0.set(2, 2, 6.0);
    try std.testing.expectApproxEqAbs(s0.get(2, 2), 6.0, std.math.floatEps(f32));
}
test "SymmetricMatrix bug test" {
    const S3f32 = SymmetricMatrix(3, f32);
    var s0 = S3f32.init(.{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
    // Setting and getting across the diagonal
    s0.set(0, 1, 99.0);
    try std.testing.expectApproxEqAbs(s0.get(1, 0), 99.0, std.math.floatEps(f32));
}

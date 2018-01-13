#include "CryptGMM/Matrix.hpp"
#include "CryptGMM/DoublePacking.hpp"
#include "HElib/EncryptedArray.h"
namespace internal 
{
static inline long round_div(long a, long b) {
    return (a + b - 1) / b;
}

/// Partition one block of the matrix into each slot of CRT packing. 
/// That is one row of the block matrix as one slot (i.e., poly).
PackedRows partition(Matrix const& matrix, BlockId const& blk, 
                     const EncryptedArray *ea, const bool backward) {
    const long d = ea->getDegree(); /// degree of the slot
    const long l = ea->size(); /// number of slots
    const long MAX_X = round_div(matrix.NumRows(), l);
    const long MAX_Y = round_div(matrix.NumCols(), d);
    assert(blk.x >= 0 && blk.x < MAX_X);
    assert(blk.y >= 0 && blk.y < MAX_Y);

    const long row_start = blk.x * l;
    const long row_end = std::min(row_start + l, matrix.NumRows());
    const long col_start = blk.y * d;
    const long col_end = std::min(col_start + d, matrix.NumCols());

    PackedRows ret;
    ret.polys.resize(l);
    ret.num_duplication = 1;
    for (long row = row_start; row < row_end; row++) {
        const long offset = row - row_start;
        ret.polys[offset].SetLength(d + 1);
        for (long col = col_start; col < col_end; col++) {
            /// coeff = backward ? d - 1 - (col - col_start) : col - col_start;
            long coeff = col - col_start;
            if (backward)
                coeff = d - 1 - coeff;
            assert(coeff >= 0 && coeff <= d - 1);
            NTL::SetCoeff(ret.polys[offset], coeff, matrix[row][col]);
        }
    }
    return ret;
}
} // namespace internal

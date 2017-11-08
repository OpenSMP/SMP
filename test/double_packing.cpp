#include <gtest/gtest.h>
#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>

#include "CryptGMM/DoublePacking.hpp"
#include "CryptGMM/Matrix.hpp"
namespace testing
{
 namespace internal
 {
    enum GTestColor {
             COLOR_DEFAULT,
             COLOR_RED,
             COLOR_GREEN,
             COLOR_YELLOW
         };
  
    extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
   }
}

#define PRINTF(...)  \
        do { \
                testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); \
                testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__); } \
                    while(0)


namespace {
const long M = 4096;
FHEcontext context(M, 113, 4);
FHESecKey *key;
Matrix mat;
class DoublePackingTest : public ::testing::Test {
protected:
    static void SetUpTestCase() {
        buildModChain(context, 6);
        key = new FHESecKey(context);
        key->GenSecKey(64);
        
        mat.SetDims(16, 16);
        for (long i = 0; i < mat.NumRows(); i++)
            for (long j = 0; j < mat.NumCols(); j++)
                mat[i][j] = j;
    }
    static void TearDownTestCase() {
        delete key;
    }
};

inline long round_div(long a, long b) {
    return (a + b - 1) / b;
}

TEST_F(DoublePackingTest, PackForward) {
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    const long MAX_X = round_div(mat.NumRows(), l);
    const long MAX_Y = round_div(mat.NumCols(), d);
    PRINTF("#slots = %ld, deg = %ld\n", l, d);
    internal::BlockId blk;
    for (int x = 0; x < MAX_X; x++) {
        for (int y = 0; y < MAX_Y; y++) {
            blk = {x, y};
            auto polys = internal::partition(mat, blk, *ea, false);

            const long row_start = blk.x * l;
            const long row_end = std::min(row_start + l, mat.NumRows());
            const long col_start = blk.y * d;
            const long col_end = std::min(col_start + d, mat.NumCols());

            for (long i = 0; i < polys.size(); i++) {
                long row = row_start + i;
                ASSERT_TRUE(row < row_end);
                for (long dd = 0; dd < d; dd++) {
                    long col = col_start + dd;
                    if (col >= col_end) 
                        break;
                    EXPECT_EQ(polys[i][dd], mat[row][col]);
                }
            }
        }
    }
}    

TEST_F(DoublePackingTest, PackBackward) {
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    const long MAX_X = round_div(mat.NumRows(), l);
    const long MAX_Y = round_div(mat.NumCols(), d);
    PRINTF("#slots = %ld, deg = %ld\n", l, d);
    internal::BlockId blk;
    for (int x = 0; x < MAX_X; x++) {
        for (int y = 0; y < MAX_Y; y++) {
            blk = {x, y};
            auto polys = internal::partition(mat, blk, *ea, true);

            const long row_start = blk.x * l;
            const long row_end = std::min(row_start + l, mat.NumRows());
            const long col_start = blk.y * d;
            const long col_end = std::min(col_start + d, mat.NumCols());

            for (long i = 0; i < polys.size(); i++) {
                long row = row_start + i;
                ASSERT_TRUE(row < row_end);
                for (long dd = 0; dd < d; dd++) {
                    long col = col_start + dd;
                    if (col >= col_end) 
                        break;
                    EXPECT_EQ(polys[i][d - dd - 1], mat[row][col]);
                }
            }
        }
    }
} 

TEST_F(DoublePackingTest, GMM) {
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    Matrix matT;
    transpose(&matT, mat);
    const long MAX_X1 = round_div(mat.NumRows(), l);
    const long MAX_Y1 = round_div(mat.NumCols(), d);
    const long MAX_X2 = MAX_Y1;
    const long MAX_Y2 = round_div(matT.NumCols(), d);
    const long plaintext = context.alMod.getPPowR();

    Matrix result = mul(mat, mat);
    for (int x = 0; x < MAX_X1; x++) {
        for (int y = 0; y < MAX_Y2; y++) {
            NewPlaintextArray sum(*ea);
            for (int k = 0; k < MAX_X2; k++) {
                internal::BlockId blk1 = {x, k};
                internal::BlockId blk2 = {k, y};
                std::vector<NTL::ZZX> polys = internal::partition(mat, blk1, *ea, true);
                NewPlaintextArray pA1(*ea); 
                encode(*ea, pA1, polys);

                polys = internal::partition(matT, blk2, *ea, false);
                NewPlaintextArray pA2(*ea);
                encode(*ea, pA2, polys);

                mul(*ea, pA1, pA2);
                add(*ea, sum, pA1);
            }

            std::vector<NTL::ZZX> decoded;
            decode(*ea, decoded, sum); 
            const long start_row = x * l;
            const long end_row = std::min(start_row + l, mat.NumRows());
            const long start_col = y * d;
            const long end_col = std::min(end_col + d, mat.NumCols());
            for (long i = 0; i < l; i++) {
                long ii = start_row + i;
                if (ii >= end_row) break;
                if (ii >= end_col) break;
                long computed = NTL::to_long(NTL::coeff(decoded[i], d - 1));
                long ground_truth = result[ii][ii] % plaintext;
                PRINTF("compute %ld, ground_truth %ld\n", computed, ground_truth);
                //EXPECT_EQ(computed, ground_truth);
            }
        }
        PRINTF("------\n");
    }
}
}

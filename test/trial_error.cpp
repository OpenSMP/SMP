#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>

#include "CryptGMM/DoublePacking.hpp"
#include "CryptGMM/Matrix.hpp"

inline long round_div(long a, long b) {
    return (a + b - 1) / b;
}

void randomize(Matrix &mat) {
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = i * 10 + j;
}

void print_diag(Matrix const& mat, int x, int height, int y, int width) {
    int x_end = std::min<int>(x + height, mat.NumRows());
    int y_end = std::min<int>(y + width, mat.NumCols());
    int window = std::min(x_end - x, y_end - y);
    printf("%d:%d;%d:%d\n", x, x_end, y, y_end);
    for (long i = 0; i < window; i++)
        std::cout << mat[x + i][y + i] << ",";
    std::cout << std::endl;
}

int main() {
    const long M = 512;
    FHEcontext context(M, 769, 4);
    FHESecKey *key;
    buildModChain(context, 6);
    key = new FHESecKey(context);
    key->GenSecKey(64);

    Matrix A, B;
    A.SetDims(200, 33);
    B.SetDims(33, 16);
    randomize(A);
    randomize(B);
    auto ground_truth = mul(A, B);

    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    std::cout << l << "," << d << std::endl;

    const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);
    /// Use the transpose of B for the sake of easy implementation.
    const long MAX_X2 = round_div(B.NumRows(), d);
    const long MAX_Y2 = round_div(B.NumCols(), l);
    const long plaintext = context.alMod.getPPowR();

    /// Use the transpose of B for the sake of easy implementation.
    Matrix Bt;
    transpose(&Bt, B);
    for (int x = 0; x < MAX_X1; x++) {
        for (int y = 0; y < MAX_Y2; y++) {
            NewPlaintextArray sum(*ea);
            for (int k = 0; k < MAX_Y1; k++) {
                internal::BlockId blk1 = {x, k};
                /// Since we use B transpose, and thus the block id is
                /// (y, k) instead of (k, y)
                internal::BlockId blk2 = {y, k};
                auto polys = internal::partition(A, blk1, *ea, false);
                NewPlaintextArray pA1(*ea); 
                encode(*ea, pA1, polys);

                polys = internal::partition(Bt, blk2, *ea, true);
                NewPlaintextArray pA2(*ea);
                encode(*ea, pA2, polys);
                mul(*ea, pA1, pA2);
                add(*ea, sum, pA1);
            }
            std::vector<NTL::ZZX> decoded;
            decode(*ea, decoded, sum); 
            for (auto& __poly : decoded)
                std::cout << NTL::to_long(NTL::coeff(__poly, d - 1)) <<",";
            std::cout << "\n-----\n";
            const long row_start = x * l;
            const long col_start = y * l;
            print_diag(ground_truth, row_start, l, col_start, l);
        }
    }
    delete key;
    return 0;
}

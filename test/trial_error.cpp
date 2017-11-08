#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>

#include "CryptGMM/DoublePacking.hpp"
#include "CryptGMM/Matrix.hpp"
#include "CryptGMM/Timer.hpp"

inline long round_div(long a, long b) {
    return (a + b - 1) / b;
}

void zero(Matrix &mat) {
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = 0;
}

void randomize(Matrix &mat) {
    for (long i = 0; i < mat.NumRows(); i++)
        for (long j = 0; j < mat.NumCols(); j++)
            mat[i][j] = NTL::RandomBnd(10L);//i * 10 + j;
}

void fill_compute(Matrix& mat, int x, int y, 
                  std::vector<NewPlaintextArray> const& inner_products,
                  const EncryptedArray *ea) {
    const long l = ea->size();
    const long d = ea->getDegree();
    /// For now only consider matrices with exact multiple l-rows and l-columns
    assert(mat.NumRows() % l == 0);
    assert(mat.NumCols() % l == 0);
    long row_start = x * l;
    long row_end = row_start + l;
    long col_start = y * l;
    long col_end = col_start + l;
    assert(row_end < mat.NumRows());
    assert(col_end < mat.NumCols());

    for (size_t i = 0; i < inner_products.size(); i++) {
        std::vector<NTL::ZZX> decoded;
        decode(*ea, decoded, inner_products[i]);
        assert(decoded.size() == l);
        for (long ll = 0; ll < l; ll++) {
            long computed = NTL::to_long(NTL::coeff(decoded[ll], d - 1));
            long row = row_start + ll;
            long col = col_start + ll + i;
            if (col >= col_end)
                col = col_start + col % l;
            assert(row < row_end);
            assert(col < col_end);
            mat.put(row, col, computed);
        }
    }
}

long ceil_round(long a, long l) {
    return (a + l - 1) / l * l;
}

int main() {
    const long M = 8192;
    //FHEcontext context(M, 769, 4);
    FHEcontext context(M, 641, 4);
    FHESecKey *key;
    buildModChain(context, 6);
    key = new FHESecKey(context);
    key->GenSecKey(64);

    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    printf("#slots = %ld, deg = %ld\n", l, d);

    Matrix A, B;
    A.SetDims(ceil_round(64, l), ceil_round(576, d));
    B.SetDims(ceil_round(576, d), ceil_round(64, l));
    randomize(A);
    randomize(B);
    printf("|A| = %ldx%ld, |B| = %ldx%ld\n", 
           A.NumRows(), A.NumCols(), B.NumRows(), B.NumCols());
    auto ground_truth = mul(A, B);

    const long MAX_X1 = round_div(A.NumRows(), l);
    const long MAX_Y1 = round_div(A.NumCols(), d);
    /// Use the transpose of B for the sake of easy implementation.
    Matrix Bt;
    transpose(&Bt, B);
    const long MAX_X2 = round_div(Bt.NumRows(), l);
    const long MAX_Y2 = round_div(Bt.NumCols(), d);
    const long plaintext = context.alMod.getPPowR();
    Matrix computed;
    computed.SetDims(A.NumRows(), B.NumCols());
    zero(computed);
    auto start_time = Clock::now(); 
    std::vector<std::vector<NewPlaintextArray>> uploading;
    uploading.resize(MAX_X1, std::vector<NewPlaintextArray>(MAX_Y1, *ea));
    for (int x = 0; x < MAX_X1; x++) {
        for (int k = 0; k < MAX_Y1; k++) {
            internal::BlockId blk = {x, k};
            auto packed_rows = internal::partition(A, blk, *ea, false);
            encode(*ea, uploading[x][k], packed_rows.polys);
        }
    }
    double encode_time = time_as_millsecond(Clock::now() - start_time);

    start_time = Clock::now();
    std::vector<std::vector<NewPlaintextArray>> precomputed;
    precomputed.resize(MAX_X2, std::vector<NewPlaintextArray>(MAX_Y1, *ea));
    for (int y = 0; y < MAX_X2; y++) {
        for (int k = 0; k < MAX_Y1; k++) {
            internal::BlockId blk = {y, k};
            auto packed_rows = internal::partition(Bt, blk, *ea, true);
            encode(*ea, precomputed[y][k], packed_rows.polys);
        }
    }
    double precompute_time = time_as_millsecond(Clock::now() - start_time);

    size_t download_ctxt = 0;
    for (int x = 0; x < MAX_X1; x++) {
        for (int y = 0; y < MAX_X2; y++) {
            std::vector<NewPlaintextArray> summations(l, *ea);
            for (int k = 0; k < MAX_Y1; k++) {
                for (long ll = 0; ll < l; ll++) {
                    auto rotated(precomputed[y][k]);
                    rotate(*ea, rotated, -ll);
                    mul(*ea, rotated, uploading[x][k]);
                    add(*ea, summations[ll], rotated);
                }
            }
            download_ctxt += summations.size();
            fill_compute(computed, x, y, summations, ea);
        }
    }

    size_t upload_ctxt = uploading.size() * uploading[0].size();
    printf("upload cipher %zd, download cipher %zd\n", upload_ctxt, download_ctxt);
    printf("%f %f\n", encode_time, precompute_time);
    if (computed != ground_truth)
        std::cerr << "Computation seems wrong" << std::endl;
    else
        std::cout << "pass" << std::endl;
    delete key;
    return 0;
}

/// Implement some functions that absent in the origin HElib.
#ifndef CRYPT_GMM_HELIB_HPP
#define CRYPT_GMM_HELIB_HPP
#include <vector>
#include <NTL/lzz_p.h>
#include <NTL/ZZX.h>
#include <memory>
class FHEcontext;
class FHESecKey;
class Ctxt;
namespace NTL { class zz_pX; }
/// encode and decode without remapping to G(X) as the EncrypedArray does.
void rawEncode(NTL::zz_pX &out, 
               std::vector<NTL::zz_pX> const& slots, 
               FHEcontext const& context);
void rawEncode(NTL::ZZX &out, 
               std::vector<NTL::zz_pX> const& slots, 
               FHEcontext const& context);
void rawDecode(std::vector<NTL::zz_pX> &out, 
               NTL::zz_pX const& poly, 
               FHEcontext const& context);
void rawDecode(std::vector<NTL::zz_pX> &out, 
               NTL::ZZX const& poly, 
               FHEcontext const& context);
void rawDecode(std::vector<NTL::ZZX> &out, 
               NTL::ZZX const& poly, 
               FHEcontext const& context);

struct GMMPrecompTable {
    std::vector<long> beta_powers; // beta^2, beta^3, ..., beta^\ell
    NTL::mulmod_t inv_p;
};

std::vector<GMMPrecompTable> precompute_gmm_tables(FHEcontext const& context);
/// Extract inner products from the decrypted polynomial
void extract_inner_products(std::vector<long> &out,
                            NTL::ZZX const& poly,
                            std::vector<GMMPrecompTable> const& tables,
                            FHEcontext const& context);

void extract_inner_products(std::vector<long> &out,
                            NTL::Vec<long> const& poly,
                            std::vector<GMMPrecompTable> const& tables,
                            FHEcontext const& context);
// This aux struct is used to extract the d-th coefficients of the packed polynomials.
struct CoeffExtractorAux {
    NTL::ZZX alpha; // from EncrypedArray::buildLinPolyCoeffs(...)
    long m; // m-th cyclotomoic poly
    long t; // plaintext Z_t
    long d; // degree of each plaintext slot
    std::vector<NTL::ZZX> merge_offsets; // packed of X^1, X^2, X^3, ..., X^{d-1}
};

bool init_coeff_extractor_aux(CoeffExtractorAux *aux, FHEcontext const& context);

// Ctxt encrypt l polynomials in the plaintext slot [A1(X), A2(X), ..., Al(X)]
// Extract the last coefficent of each packed polynomial in-place.
void extract_last_coeffient(Ctxt &ctxt, CoeffExtractorAux const& aux);
// merge cnt <= d ciphertexts as a single one, by shifting the coefficients and then summing up the shifted ciphertexts.
// i.e., ctx(X) + ctx(X^2) + ctx(X^3) + G
Ctxt merge_ctxts_by_shifting(Ctxt const* ctxt, size_t cnt, CoeffExtractorAux const& aux);
#endif //CRYPT_GMM_HELIB_HPP

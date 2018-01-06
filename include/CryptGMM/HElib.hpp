/// Implement some functions that absent in the origin HElib.
#ifndef CRYPT_GMM_HELIB_HPP
#define CRYPT_GMM_HELIB_HPP
#include <vector>
#include <NTL/lzz_p.h>
class FHEcontext;
namespace NTL { class zz_pX; class ZZX; }
/// encode and decode without using the G(X) as the EncrypedArray does.
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
    std::vector<long> beta_powers;
    NTL::mulmod_t inv_p;
};

std::vector<GMMPrecompTable> precompute_gmm_tables(FHEcontext const& context);
/// Extract inner products from the decrypted polynomial
void extract_inner_products(std::vector<long> &out,
                            NTL::ZZX const& poly,
                            std::vector<GMMPrecompTable> const& tables,
                            FHEcontext const& context);

#endif //CRYPT_GMM_HELIB_HPP

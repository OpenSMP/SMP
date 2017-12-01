/// Implement some functions that absent in the origin HElib.
#ifndef CRYPT_GMM_HELIB_HPP
#define CRYPT_GMM_HELIB_HPP
#include <vector>
class FHEcontext;
namespace NTL { class zz_pX; class ZZX; }
/// encode and decode without using the G(X) as the EncrypedArray does.
void rawEncode(NTL::zz_pX &out, std::vector<NTL::zz_pX> const& slots, FHEcontext const& context);
void rawEncode(NTL::ZZX &out, std::vector<NTL::zz_pX> const& slots, FHEcontext const& context);
void rawDecode(std::vector<NTL::zz_pX> &out, NTL::zz_pX const& poly, FHEcontext const& context);
void rawDecode(std::vector<NTL::zz_pX> &out, NTL::ZZX const& poly, FHEcontext const& context);
void rawDecode(std::vector<NTL::ZZX> &out, NTL::ZZX const& poly, FHEcontext const& context);
#endif //CRYPT_GMM_HELIB_HPP

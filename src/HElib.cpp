#include "CryptGMM/HElib.hpp"
#include <HElib/FHEContext.h>
#include <HElib/PAlgebra.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
void rawEncode(NTL::zz_pX &out,
               std::vector<NTL::zz_pX> const& slots,
               FHEcontext const& context)
{
    const auto &encoder = context.alMod.getDerived(PA_zz_p());
    encoder.CRT_reconstruct(out, const_cast<std::vector<NTL::zz_pX> &>(slots));
}

void rawEncode(NTL::ZZX &out,
               std::vector<NTL::zz_pX> const& slots,
               FHEcontext const& context)
{
    NTL::zz_pX tmp;
    rawEncode(tmp, slots, context);
    NTL::conv(out, tmp);
}

void rawDecode(std::vector<NTL::zz_pX> &out,
               NTL::zz_pX const& poly,
               FHEcontext const& context)
{
    const auto &decoder = context.alMod.getDerived(PA_zz_p());
    decoder.CRT_decompose(out, poly);
}

void rawDecode(std::vector<NTL::zz_pX> &out,
               NTL::ZZX const& poly,
               FHEcontext const& context)
{
    NTL::zz_pX tmp;
    NTL::conv(tmp, poly);
    rawDecode(out, tmp, context);
}

void rawDecode(std::vector<NTL::ZZX> &out,
               NTL::ZZX const& poly,
               FHEcontext const& context)
{
    NTL::zz_pX tmp;
    std::vector<NTL::zz_pX> tmp2;
    rawDecode(tmp2, tmp, context);
    out.resize(tmp2.size());
    for (size_t i = 0; i < tmp2.size(); i++)
        NTL::conv(out[i], tmp2[i]);
}

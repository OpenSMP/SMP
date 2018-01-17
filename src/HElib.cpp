#include "SMP/HElib.hpp"
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/PAlgebra.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>

#include <algorithm>
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

static GMMPrecompTable precompute_gmm_table(long beta, long p, long slots)
{
    GMMPrecompTable tbl;
    tbl.beta_powers.resize(slots);
    /// (-beta)^k mod p for 0 <= k < l
    for (long i = 0; i < slots; i++)
        tbl.beta_powers[i] = NTL::PowerMod(i & 1 ? p - beta : beta, i, p);
    return tbl;
}

static bool is_valid_for_GMM(NTL::ZZX const& factor)
{
    for (long d = 1; d < NTL::deg(factor); d++) {
        if (NTL::coeff(factor, d) != 0)
            return false;
    }
    return true;
}

std::vector<GMMPrecompTable> precompute_gmm_tables(FHEcontext const& context)
{
    long p = context.alMod.getPPowR();
    long l = context.ea->size();
    std::vector<GMMPrecompTable> tbls;
    tbls.reserve(l);
    NTL::mulmod_t inv_p = NTL::PrepMulMod(p);
    for (const auto & factor : context.alMod.getFactorsOverZZ()) {
        assert(is_valid_for_GMM(factor));
        long beta = NTL::to_long(factor[0]);
        tbls.push_back(precompute_gmm_table(beta, p, l)) ;
        tbls.back().inv_p = inv_p;
    }
    return tbls;
}

long extract_inner_product(NTL::ZZX const& poly,
                           GMMPrecompTable const& tbl,
                           FHEcontext const& context)
{
    long d = context.ea->getDegree();
    long l = context.ea->size();
    long p = context.alMod.getPPowR();
    long phim = context.zMStar.getPhiM();
    long ret = 0;
    for (long i = 0; i < l; i++) {
        long coeff_loc = (i + 1) * d - 1;
        assert(coeff_loc < phim);
        long coeff = NTL::to_long(NTL::coeff(poly, coeff_loc));
        coeff = NTL::MulMod(coeff, tbl.beta_powers.at(i), p, tbl.inv_p);
        ret = NTL::AddMod(ret, coeff, p);
    }
    return ret;
}

long extract_inner_product(NTL::Vec<long> const& poly,
                           GMMPrecompTable const& tbl,
                           FHEcontext const& context)
{
    long d = context.ea->getDegree();
    long l = context.ea->size();
    long p = context.alMod.getPPowR();
    long phim = context.zMStar.getPhiM();
    assert(poly.length() == phim);
    long ret = 0;
    for (long i = 0; i < l; i++) {
        long coeff_loc = (i + 1) * d - 1;
        assert(coeff_loc < phim);
        long coeff = poly.at(coeff_loc) % p;
        coeff = NTL::MulMod(coeff, tbl.beta_powers.at(i), p, tbl.inv_p);
        ret = NTL::AddMod(ret, coeff, p);
    }
    return ret;
}

void extract_inner_products(std::vector<long> &out,
                            NTL::ZZX const& poly,
                            std::vector<GMMPrecompTable> const& tables,
                            FHEcontext const& context)
{
    out.clear();
    out.reserve(context.ea->size());
    for (const auto &tbl : tables) {
        out.push_back(extract_inner_product(poly, tbl, context));
    }
}

void extract_inner_products(std::vector<long> &out,
                            NTL::Vec<long> const& poly,
                            std::vector<GMMPrecompTable> const& tables,
                            FHEcontext const& context)
{
    out.clear();
    out.reserve(context.ea->size());
    for (const auto &tbl : tables) {
        out.push_back(extract_inner_product(poly, tbl, context));
    }
}

void faster_decrypt(NTL::Vec<long> &out,
                    FHESecKey const& sk,
                    Ctxt const& ctx)
{
    const FHEcontext& context = sk.getContext();
    assert(ctx.getPrimeSet().card() == 1);
    DoubleCRT dcrt(context, ctx.getPrimeSet()); // Set to zero
    sk.Decrypt(dcrt, ctx);
    dcrt.getOneRow(out, 0);
}

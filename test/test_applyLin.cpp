//
// Created by riku on 6/20/18.
//
#include "HElib/FHE.h"
#include "HElib/FHEContext.h"
#include "HElib/PAlgebra.h"
#include "HElib/EncryptedArray.h"
#include <NTL/ZZX.h>

#include "SMP/HElib.hpp"
#include "SMP/Timer.hpp"

void frobenious(
    NTL::zz_pX &H,
    long j,
    const FHEcontext &context)
{
    long t = context.zMStar.getP();
    long d = context.ea->getDegree();
    long jj = mcMod(j, d);
    NTL::zz_pX G;
    NTL::conv(G, context.alMod.getFactorsOverZZ()[0]);
    auto tmp = PowerMod(NTL::zz_pX(1, 1), NTL::power_ZZ(t, jj), G);
    H = CompMod(H, tmp, G);
}

void MyBuildLinPoly(std::vector<std::vector<NTL::ZZX>> &Cs,
                    const std::vector<NTL::ZZX> &Ls,
                    const EncryptedArray *ea)
{
    using RBak = PA_zz_p::RBak;
    using REbak = PA_zz_p::REBak;
    using RE = PA_zz_p::RE;
    using RX = PA_zz_p::RX;

    RBak bak; bak.save(); ea->restoreContext();
    REbak ebak; ebak.save(); ea->restoreContextForG();
    NTL::Lazy< NTL::Mat<RE> > linPolyMatrix;
    do {
        typename Lazy< Mat<RE> >::Builder builder(linPolyMatrix);
        if (!builder()) break;

        FHE_NTIMER_START(buildLinPolyCoeffs_invert);


        long p = ea->getPAlgebra().getP();
        long r = ea->getAlMod().getR();

        NTL::Mat<RE> M1;
        // build d x d matrix, d is taken from the current NTL context for G
        buildLinPolyMatrix(M1, p);
        NTL::Mat<RE> M2;
        ppInvert(M2, M1, p, r); // invert modulo prime-power p^r

        NTL::UniquePtr< Mat<RE> > ptr;
        ptr.make(M2);
        builder.move(ptr);
    } while (0);

    long d = ea->getDegree();
    Cs.resize(d);
    for (long i = 0; i < d; ++i) {
        RE L = NTL::conv<RE>(NTL::conv<RX>(Ls.at(i)));
        NTL::Vec<RE> tmp;
        NTL::mul(tmp, L, linPolyMatrix->operator[](d - 1));
        convert(Cs.at(i), tmp);
    }
}

#if 0
int main2() {
    long m = 8192 << 1;
    long t = 67073;
    FHEcontext context(m, t, 1);
    context.bitsPerLevel = 60;
    buildModChain(context, 4);
    const EncryptedArray *ea = context.ea;
    std::cout << "G = " << context.alMod.getFactorsOverZZ()[0] << std::endl;

    std::vector<NTL::ZZX> slots(ea->size());
    for (auto &&ss : slots) {
        ss.SetLength(ea->getDegree());
        for (long d = 0; d < ea->getDegree(); ++d)
            NTL::SetCoeff(ss, d, d + 1);
    }

    const long l = ea->size();
    const long d = ea->getDegree();
    std::vector<NTL::ZZX> L(ea->getDegree(), NTL::ZZX(0));
    L.back() = NTL::ZZX(0, 1);

    std::vector<NTL::ZZX> alphas;
    ea->buildLinPolyCoeffs(alphas, L);
    std::cout << alphas << "\n";

    CoeffExtractorAux aux;
    {
        NTL::ZZX alpha;
        MyBuildLinPoly(alpha,L.back(), ea);
        std::vector<NTL::ZZX> tmp(ea->size(), alpha);
        ea->encode(aux.alpha, tmp);
    }
    aux.m = m;
    aux.t = t;
    aux.d = ea->getDegree();

    FHESecKey sk(context);
    sk.GenSecKey(64);
    addFrbMatrices(sk);

    Ctxt ctx(sk);
    NTL::ZZX plain;
    ea->encode(plain, slots);
    sk.Encrypt(ctx, plain);
    ctx.multByConstant(plain);

    Ctxt cpy(ctx);
    double extract_time = 0.;
    {
        AutoTimer timer(&extract_time);
        extract_last_coeffient(ctx, aux);
    }
    std::cout << "time "  << extract_time << "\n";

    extract_time = 0.;
    {
        AutoTimer timer(&extract_time);
        applyLinPoly1(*ea, cpy, alphas);
    }
    std::cout << "time "  << extract_time << "\n";
    std::cout << ctx.log_of_ratio() << ":" << cpy.log_of_ratio() << std::endl;

    ea->decrypt(ctx, sk, slots);
    std::cout << slots.back() << "\n";

    ea->decrypt(cpy, sk, slots);
    std::cout << slots.back() << "\n";
    return 0;
}
#endif

void create_coeff(
     std::vector<NTL::ZZX> &coeff,
     long pos,
     const EncryptedArray *ea)
{
    std::vector<NTL::ZZX> L(ea->getDegree(), NTL::ZZX(0));
    L.back() = NTL::ZZX(pos, 1);
    ea->buildLinPolyCoeffs(coeff, L);
}

template <class T>
double obj_size(T const& obj)
{
    std::stringstream ss;
    ss << obj;
    std::string s = ss.str();
    return s.size() / 1024.;
}

int main() {
    //long m = 8192;
    //long t = 67073;
    long m = 16384;
    long t = 70657;
    FHEcontext context(m, t, 1);
    context.bitsPerLevel = 50;
    buildModChain(context, 3);
    const EncryptedArray *ea = context.ea;
    FHESecKey sk(context);
    sk.GenSecKey(64, 0, 1);
    addFrbMatrices(sk);
    FHEPubKey const &pk{sk};
    const long l = ea->size();
    const long d = ea->getDegree();

    NTL::zz_p::init(t);
    std::vector<zzX> slots(l);
    for (auto &&ss : slots) {
        ss.SetLength(d);
        for (long i = 0; i < d; ++i)
            ss[i] = i + 1;
    }

    std::vector<std::vector<NTL::ZZX>> coeffs(d);
    std::vector<NTL::ZZX> Ls(d, NTL::ZZX(0));
    double prepare_time = 0.;
    {
        AutoTimer timer(&prepare_time);
        for (long i = 0; i < d; ++i)
            Ls[i] = NTL::ZZX(i, 1);
        MyBuildLinPoly(coeffs, Ls, ea);

        for (long i = 0; i < d; ++i) {
            size_t sze = coeffs.at(i).size();
            std::vector<NTL::ZZX> encoded_coeffs(sze);
            for (long j = 0; j < sze; ++j) {
                std::vector<NTL::ZZX> tmp(l, coeffs[i][j]);
                ea->encode(encoded_coeffs[j], tmp);
            }
            std::swap(coeffs[i], encoded_coeffs);
        }
    }
    std::cout << "prepare all linear maps took " << prepare_time << " ms\n";

    Ctxt ctx(pk);
    zzX _plain;
    ea->encode(_plain, slots);
    NTL::ZZX plain;
    convert(plain, _plain);
    sk.Encrypt(ctx, plain);
    ctx.multByConstant(_plain);
    for (long i = 0; i < 10; ++i)
        ctx += ctx;
    Ctxt origin(ctx);

    std::vector<Ctxt> ctxs(d, ctx);
    double extra_time = 0;
    double merge_time = 0.;
    {
        AutoTimer timer(&extra_time);
        for (long i = 0; i < d; ++i)
            applyLinPolyLL(ctxs[i], coeffs[i], d);
    }
    {
        AutoTimer timer(&merge_time);
        for (long i = 1; i < d; ++i)
            ctxs[0] += ctxs[i];
    }
    ctxs[0].modDownToLevel(1);

    double dec_time = 0.;
    std::vector<NTL::ZZX> decoded;
    {
        AutoTimer timer(&dec_time);
        NTL::ZZX dec;
        sk.Decrypt(dec, ctxs[0]);
        ea->decode(decoded, dec);
    }
    std::cout << decoded << std::endl;

    std::cout << "security level " << context.securityLevel() << std::endl;
    std::cout << "d@l = " << d << "@" << l << std::endl;
    std::cout << "|evk| = " << obj_size(sk) << std::endl;

    std::cout << "final noise " << ctxs[0].log_of_ratio() << std::endl;
    std::cout << "|ctx_0| = " << obj_size(origin) << std::endl;
    std::cout << "|ctx_L| = " << obj_size(ctxs[0]) << std::endl;

    std::cout << "extract " << d << " ctxts took " << extra_time << "ms" << std::endl;
    std::cout << "merge " << d << " ctxts took " << merge_time << "ms" << std::endl;
    std::cout << "one dec took " << dec_time << "ms" << std::endl;

    return 0;
}

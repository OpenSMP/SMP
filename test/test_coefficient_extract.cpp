#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/PAlgebra.h>
#include <NTL/ZZX.h>

#include "SMP/HElib.hpp"
#include "SMP/Timer.hpp"

Ctxt merge(std::vector<Ctxt> arr, const EncryptedArray *ea)
{
    assert(!arr.empty());
    size_t sze = arr.size();
    const long l = ea->size();
    Ctxt ans{arr.front()};
    for (long j = 1; j < sze; ++j) {
        std::vector<NTL::ZZX> slots(l, NTL::ZZX(j, 1));
        NTL::ZZX offset;
        ea->encode(offset, slots);
        arr[j].multByConstant(offset);
        ans += arr[j];
    }
    return ans;
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
    long m = 8192 << 1;
    long t = 67073;
    FHEcontext context(m, t, 1);
    context.bitsPerLevel = 60;
    buildModChain(context, 4);
    std::cout << "|g| = " << context.zMStar.numOfGens() << "\n";
    std::cout << "kappa = " << context.securityLevel() << "\n";
    FHESecKey sk(context);
    sk.GenSecKey(64, 0, 1); // no need s^2 -> s KW
    addFrbMatrices(sk);
    const FHEPubKey &pk{sk};
    std::cout << "|evk| = " << obj_size(pk) << "\n";
    
    const EncryptedArray *ea = context.ea;
    const long d = ea->getDegree();
    const long l = ea->size();
    std::cout << "d@l = " << d << "@" << l << std::endl;

    NTL::ZZX poly;
    //NTL::zz_pX poly;
    poly.SetLength(d + 1);
    for (long i = 0; i < d; ++i)
        NTL::SetCoeff(poly, i, i + 1);

    std::vector<NTL::ZZX> slots(l, poly);
    NTL::ZZX plain;
    ea->encode(plain, slots);

    Ctxt ctx(sk);
    sk.Encrypt(ctx, plain);
    std::cout << "|ctxt| = " << obj_size(ctx) << "\n";
    double multByConstant = 0.;
    {
        AutoTimer timer(&multByConstant);
        ctx.multByConstant(plain); // mimic the vector * vector
    }

    CoeffExtractorAux aux;
    init_coeff_extractor_aux(&aux, context);
    double extract = 0.;
    {
        AutoTimer timer(&extract);
        extract_last_coeffient(ctx, aux); // extract the inner products
    }

    std::vector<Ctxt> copies(d, ctx); // mimic a collection of inner products
    double merge_time = 0.;
    {
        AutoTimer timer(&merge_time);
        ctx = merge_ctxts_by_shifting(copies.data(), copies.size(),
                                      aux); // merge multiple inner products in a single ctxt
        ctx.modDownToLevel(1);
    }
    std::cout << "|ctxt| = " << obj_size(ctx) << "\n";

    std::cout << "Correct ? " << ctx.isCorrect() << std::endl;
    std::cout << "log ratio: " << ctx.log_of_ratio() << "\n";
    double dec_time = 0., decode_time = 0.;
    {
        AutoTimer timer(&dec_time);
        sk.Decrypt(plain, ctx);
    }

    {
        AutoTimer timer(&decode_time);
        ea->decode(slots, plain);
    }

    std::cout << slots.front() << std::endl;
    std::cout << slots.back() << std::endl;
    std::cout << multByConstant << "," <<
                 extract << "," <<
                 merge_time << "," <<
                 dec_time << "," <<
                 decode_time << "," <<
                 std::endl;
    printAllTimers(std::cout);
    return 0;
}

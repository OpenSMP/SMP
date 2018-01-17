#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/PAlgebra.h>
#include <NTL/ZZX.h>
#include "SMP/Timer.hpp"
#include "SMP/literal.hpp"
#include "SMP/HElib.hpp"
bool check(std::vector<ZZX> const& factors, long d) {
    std::cout << "factors: \n";
    for (auto const &f : factors) {
        std::cout << f << "\n";
        for (long i = 1; i < d; i++) {
            if (NTL::coeff(f, i) != 0)
                return false;
        }
    }
    return true;
}

std::ostream& print_poly(NTL::ZZX const& poly, long d) {
    std::cout << "[";
    for (; d > 0; d--)
        std::cout << NTL::coeff(poly, d) << ",";
    std::cout << NTL::coeff(poly, 0) << "]";
    return std::cout;
}

long inner_product(NTL::zz_pX const& a, 
                   NTL::zz_pX const& b) 
{
    NTL::zz_p ip;
    long deg = NTL::deg(a);
    for (long i = 0; i <= deg; i++) {
        ip += (NTL::coeff(a, i) * NTL::coeff(b, deg - i));
    }
    return ip._zz_p__rep;
}

/// Return the leading coeff of the rem polynomial f % h
/// h = X^d + a
NTL::zz_p ModGetLeadingCoeff(NTL::zz_pX const& f,
                             NTL::zz_pX const& h) 
{
    long m = NTL::deg(f) + 1;
    long d = NTL::deg(h);
    long n = m / d;

    // printf("|f| = %ld, |h| = %ld\n", m, d);
    // std::cout << f << "\n" << h << std::endl;
    std::cout << h << std::endl;
    NTL::zz_p beta = NTL::coeff(h, 0);
    NTL::zz_p beta_power(1);
    NTL::zz_p ret(NTL::coeff(f, d - 1));
    std::cout << ret;
    for (long i = 2; i <= n; i++) {
        auto coeff = NTL::coeff(f, i * d - 1);
        beta_power *= beta;
        coeff *= beta_power;
        if ((i & 1) == 0) {
            ret -= coeff;
        std::cout << " - " << beta_power << " * " << NTL::coeff(f, i * d - 1) ;
        } else {
        std::cout << " + " << beta_power << " * " << NTL::coeff(f, i * d - 1) ;
            ret += coeff;
        }
    }
    std::cout << "\n";
    auto be = beta;
    for (long i = 1; i <= n; i++) {
        std::cout << be << " ";
        be *= beta;
    }
    std::cout << "\n";
    // std::cout << "result " << ret << std::endl;
    return ret;
}

std::vector<long> precompute_power(long beta, long p, long l)
{
    std::vector<long> beta_power_l(l);
    /// (-beta)^k mod p for 0 <= k < l
    for (long i = 0; i < l; i++)
        beta_power_l[i] = NTL::PowerMod(i & 1 ? p - beta : beta, i, p);
    return beta_power_l;
}

long mod_with_precomputed_table(NTL::ZZX const& poly,
                                std::vector<long> const& tbl,
                                FHEcontext const &context)
{
    long d = context.ea->getDegree();
    long l = context.ea->size();
    long p = context.alMod.getPPowR();
    long phim = context.zMStar.getPhiM();
    long ret = 0;
    auto inv_p = NTL::PrepMulMod(p);
    for (long i = 0; i < l; i++) {
        long coeff_loc = (i + 1) * d - 1;
        assert(coeff_loc < phim);
        long coeff = NTL::to_long(NTL::coeff(poly, coeff_loc));
        coeff = NTL::MulMod(coeff, tbl[i], p, inv_p);
        ret = NTL::AddMod(ret, coeff, p);
    }
    return ret;
}

void faster() {
    NTL::SetSeed(NTL::to_ZZ(132));
    long m = 8192;
    long p = 70913;
    NTL::zz_p::init(p);
    FHEcontext context(m, p, 1);
    context.bitsPerLevel = 60;
    buildModChain(context, 2);
    FHESecKey sk(context);
    sk.GenSecKey(64);

    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    std::cout << "l = " << l << ", d = " << d << std::endl;

    std::vector<GMMPrecompTable> tables = precompute_gmm_tables(context);
    const auto &encoder = context.alMod.getDerived(PA_zz_p());

    std::vector<NTL::zz_pX> vec_A(l);
    std::vector<NTL::zz_pX> vec_B(l);
    for (long i = 0; i < l; i++) {
        NTL::random(vec_A[i], d);
        NTL::random(vec_B[i], d);
        std::cout << inner_product(vec_A[i], vec_B[i]) << " ";
    }
    std::cout << "\n";

    NTL::zz_pX encoded_A, encoded_B;
    encoder.CRT_reconstruct(encoded_A, vec_A);
    encoder.CRT_reconstruct(encoded_B, vec_B);

    Ctxt ctx_vec_A(sk), ctx_vec_B(sk);
    sk.Encrypt(ctx_vec_A, NTL::conv<NTL::ZZX>(encoded_A));
    sk.Encrypt(ctx_vec_B, NTL::conv<NTL::ZZX>(encoded_B));
    ctx_vec_A.multiplyBy(ctx_vec_B);

    std::vector<NTL::zz_pX> slots;
    NTL::ZZX decrypted;
    sk.Decrypt(decrypted, ctx_vec_A);
    double raw_dec_time = 0.;
    {
        AutoTimer timer(&raw_dec_time);
        rawDecode(slots, decrypted, context);
    }

    double table_dec_time = 0.;
    std::vector<long> inner_products;
    {
        AutoTimer timer(&table_dec_time);
        extract_inner_products(inner_products, decrypted, tables, context);
    }
    std::cout << inner_products << std::endl;

    std::cout << raw_dec_time << " : " << table_dec_time << std::endl;
}

void test_poly_mod()
{
    NTL::zz_p::init(769);
    NTL::zz_pX ply;
    ply.SetLength(3);
    ply[0] = 3; ply[2] = 1; // X^2 + 3
    NTL::zz_pXModulus mod(ply);

    NTL::zz_pX rnd;
    NTL::random(rnd, 6);
    std::cout << rnd << std::endl;

    NTL::zz_pX rm;
    NTL::rem(rm, rnd, mod);
    std::cout << rm << std::endl;
}

void normal() {
    long m = 8192;
    long p = 769;
    FHEcontext context(m, p, 1);
    context.bitsPerLevel = 59;
    buildModChain(context, 2);
    FHESecKey sk(context);
    sk.GenSecKey(64);

    const auto &factors = context.alMod.getFactorsOverZZ();
    EncryptedArray *ea = new EncryptedArray(context, factors[0]);
    const long l = ea->size();
    const long d = ea->getDegree();

    NTL::zz_p::init(p);
    std::vector<NTL::zz_pX> vec_A(l);
    std::vector<NTL::zz_pX> vec_B(l);
    std::vector<NTL::ZZX> Vec_A(l);
    std::vector<NTL::ZZX> Vec_B(l);
    for (long i = 0; i < l; i++) {
        NTL::random(vec_A[i], d);
        NTL::random(vec_B[i], d);
        NTL::conv(Vec_A[i], vec_A[i]);
        NTL::conv(Vec_B[i], vec_B[i]);
        std::cout << inner_product(vec_A[i], vec_B[i]) << " ";
    }
    std::cout << "\n";

    NTL::ZZX encoded_A, encoded_B;
    for (long i = 0; i < 100; i++) {
        FHE_NTIMER_START(EAEncode);
        ea->encode(encoded_A, Vec_A);
        ea->encode(encoded_B, Vec_B);
        FHE_NTIMER_STOP(EAEncode);

        FHE_NTIMER_START(rawEncode);
        rawEncode(encoded_A, vec_A, context);
        rawEncode(encoded_B, vec_B, context);
        FHE_NTIMER_STOP(rawEncode);
    }


    Ctxt ctx_vec_A(sk), ctx_vec_B(sk);
    sk.Encrypt(ctx_vec_A, encoded_A);
    sk.Encrypt(ctx_vec_B, encoded_B);
    ctx_vec_A.multiplyBy(ctx_vec_B);

    NTL::ZZX decrypted;
    sk.Decrypt(decrypted, ctx_vec_A);
    std::vector<NTL::zz_pX> results;
    std::vector<NTL::ZZX> Results;
    std::vector<double> raw_decodes;
    for (long i = 0; i < 100; i++) {
        FHE_NTIMER_START(EADecode);
        ea->decode(Results, decrypted);
        FHE_NTIMER_STOP(EADecode);

        FHE_NTIMER_START(rawDecode);
        do {
            raw_decodes.push_back(0.);
            AutoTimer timer(&(raw_decodes.back()));
            rawDecode(results, decrypted, context);
        } while (0);
        FHE_NTIMER_STOP(rawDecode);
    }
    //for (auto &s : results) {
    //    std::cout << NTL::coeff(s, d - 1) << " ";
    //}
    //std::cout << "\n";

    printNamedTimer(std::cout, "EAEncode");
    printNamedTimer(std::cout, "EADecode");
    printNamedTimer(std::cout, "rawEncode");
    printNamedTimer(std::cout, "rawDecode");
    auto decode = mean_std(raw_decodes);
    std::cout << decode.first / 1000. << " +- " << decode.second << std::endl;
}

int main() {
    //normal();
    faster();
    //test_poly_mod();
    return 0;
}

#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/PAlgebra.h>
#include <NTL/ZZX.h>
#include "CryptGMM/HElib.hpp"
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

void faster() {
    long m = 32;
    long p = 401;
    FHEcontext context(m, p, 1);
    buildModChain(context, 4);
    FHESecKey sk(context);
    sk.GenSecKey(64);

    std::vector<NTL::zz_pX> factors; 
    for (auto const& f: context.alMod.getFactorsOverZZ())  {
        factors.emplace_back(NTL::conv<NTL::zz_pX>(f));
    }
    const EncryptedArray *ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
    std::cout << "l = " << l << ", d = " << d << std::endl;
    const auto &encoder = context.alMod.getDerived(PA_zz_p());

    NTL::zz_p::init(p);
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

    NTL::ZZX decrypted;
    sk.Decrypt(decrypted, ctx_vec_A);
    NTL::zz_pX f = NTL::conv<NTL::zz_pX>(decrypted);
    std::cout << "f = " << f << "\n";
    std::cout << ModGetLeadingCoeff(f, factors[1]) << std::endl;
    std::cout << ModGetLeadingCoeff(f, factors[0]) << std::endl;
    // for (auto &h : factors) {
    //      std::cout << ModGetLeadingCoeff(f, h) << " " ;
    // }
    // encoder.CRT_decompose(vec_A, 
    // for (auto &s : vec_A) {
    //     std::cout << NTL::coeff(s, d - 1) << " ";
    // }
    std::cout << "\n";
}

void normal() {
    long m = 8192;
    long p = 113;
    FHEcontext context(m, p, 1);
    buildModChain(context, 4);
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
    for (long i = 0; i < 100; i++) {
        FHE_NTIMER_START(EADecode);
        ea->decode(Results, decrypted);
        FHE_NTIMER_STOP(EADecode);

        FHE_NTIMER_START(rawDecode);
        rawDecode(results, decrypted, context);
        FHE_NTIMER_STOP(rawDecode);
    }
    for (auto &s : results) {
        std::cout << NTL::coeff(s, d - 1) << " ";
    }
    std::cout << "\n";

    printNamedTimer(std::cout, "EAEncode");
    printNamedTimer(std::cout, "EADecode");
    printNamedTimer(std::cout, "rawEncode");
    printNamedTimer(std::cout, "rawDecode");
}

int main() {
    normal();
    // std::cout << "----------" << std::endl;
    // faster();
    return 0;
}

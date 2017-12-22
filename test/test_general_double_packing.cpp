#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <iostream>
#include <fstream>
void generate_poly(std::vector<NTL::ZZX> &polys,
                   const EncryptedArray *ea,
                   const long upper) {
    polys.resize(ea->size());
    long d = ea->getDegree();
    for (size_t i = 0; i < polys.size(); i++) {
        polys[i].SetLength(d);
        for (long j = 0; j < d; j++)
            NTL::SetCoeff(polys[i], j, NTL::RandomBnd(upper));
        // std::cout << NTL::deg(polys[i]) << "=?" << d << std::endl;
    }
}

long inner_prod(NTL::ZZX const& v, NTL::ZZX const& u, 
                const long d,
                const long upper) {
    NTL::ZZ ip(0);
    // std::cout << v << " " << u << std::endl;
    for (long i = 0; i < d; i++) {
        // std::cout << NTL::coeff(v, i)  << "*" << NTL::coeff(u, d - i - 1) << std::endl;
        ip = ip + NTL::coeff(v, i) * NTL::coeff(u, d - i - 1);
    }
    return ip % upper;
}

int main(int argc, char *argv[]) {
    // long m = 32;
    // long p = 113;
    long m = NTL::RandomPrime_long(5, 20);
    long p = NTL::RandomPrime_long(13, 20);
    do {
          p = NTL::RandomPrime_long(13, 20);
     } while ((m % p) == 0);
    FHEcontext context(m, p, 1);
    auto G = context.alMod.getFactorsOverZZ()[0];
    // EncryptedArray *ea = new EncryptedArray(context, G);
    auto ea = context.ea;
    std::cout << m << " " << p << std::endl;
    std::cout << ea->size() << " " << ea->getDegree() << std::endl;
    buildModChain(context, 4);
    FHESecKey sk(context);
    sk.GenSecKey(64);
    
    std::vector<NTL::ZZX> slots, slots2;
    generate_poly(slots, ea, p);
    generate_poly(slots2, ea, p);
    for (size_t i = 0; i < slots.size(); i++) {
        std::cout << inner_prod(slots[i], slots2[i], ea->getDegree(), p) << std::endl;
    }
    std::cout << "\n";
    Ctxt ctx(sk), ctx2(sk);
    ea->skEncrypt(ctx, sk, slots); 
    ea->skEncrypt(ctx2, sk, slots2); 
    ctx *= ctx2;
    ea->decrypt(ctx, sk, slots);
    for (auto &ss : slots)
        std::cout << NTL::coeff(ss, ea->getDegree() - 1) << std::endl;
    return 0;
}

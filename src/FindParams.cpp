//
// Created by riku on 1/3/18.
//
#include <NTL/ZZ.h>
#include <HElib/NumbTh.h>
#include <HElib/FHEContext.h>
#include <iostream>
void MiniONN(long m) {
    long p = NTL::RandomPrime_long(10, 20);
    while (true) {
        long d = multOrd(p, m);
        if (d == 1)
            break;
        do {
            p = NTL::NextPrime(p + 2);
        } while (m % p == 0);
    };
    std::cout << m << " " << p << std::endl;
}

bool check(NTL::ZZX const& factor) {
	for (long i = 1; i < NTL::deg(factor); i++) {
		if (NTL::coeff(factor, i) != 0)
			return false;
	}
	return true;
}

void DoublePacking(long m, long slots) {
	long p = 3;
	long phim = phi_N(m);
	assert(phim == (m >> 1));
	while (true) {
		long d = multOrd(p, m);
		long s = phim / d;
		if (s == slots) {
			FHEcontext context(m, p, 1);
			const auto &ftrs = context.alMod.getFactorsOverZZ();
			bool ok = true;
			for (const auto& f : ftrs)
				ok &= check(f);
			if (ok)
				printf("%ld %ld\n", m, p);
		}
		do {
            p = NTL::NextPrime(p + 2);
        } while (m % p == 0);
	};
}

int main() {
    //MiniONN(8192);
	DoublePacking(8192, 128);
    return 0;
}


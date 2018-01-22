#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/NumbTh.h>
#include <HElib/EncryptedArray.h>

#include "SMP/DoublePacking.hpp"
#include "SMP/Timer.hpp"
#include "SMP/HElib.hpp"

#include <iostream>
#include <vector>
void randomize(NTL::ZZX &poly, long phim, long p = 3) {
	for (long d = 0; d < phim; d++)
		NTL::SetCoeff(poly, d, NTL::RandomBnd(p));
}

double slower_version(std::vector<NTL::zz_pX> &slots,
					  NTL::ZZX const& poly, FHEcontext const& context) {
	NTL::zz_pX poly_p;
	NTL::conv(poly_p, poly);
	rawDecode(slots, poly_p, context);
}

double faster_version(std::vector<long> &slots,
					  NTL::ZZX const& poly, 
					  std::vector<GMMPrecompTable> const& tbls,
					  FHEcontext const& context) {
	extract_inner_products(slots, poly, tbls, context);
}

int run(long m, long p) {
    const long r = 1;
    const long L = 2;
    NTL::zz_p::init(p);
    FHEcontext context(m, p, r);
    context.bitsPerLevel = 60;
    buildModChain(context, L);
    std::cerr << "kappa = " << context.securityLevel() << std::endl;
    std::cerr << "slot = " << context.ea->size() << std::endl;
    std::cerr << "degree = " << context.ea->getDegree() << std::endl;

	NTL::ZZX rnd;
	rnd.SetLength(context.zMStar.getPhiM());
	randomize(rnd, (m >> 1), p);
	double pre_time = 0.;
	std::vector<GMMPrecompTable> tbls;
	{
		AutoTimer timer(&pre_time);
		tbls = precompute_gmm_tables(context);
	}

	std::vector<long> slots;
	const long T = 2;
	double faster_time = 0.;
	for (long i = 0; i < T; i++) {
		double _t = 0.;
		{
			AutoTimer timer(&_t);
			faster_version(slots, rnd, tbls, context);
		}
		faster_time += _t;
	}

	std::vector<NTL::zz_pX> _slots;
	double slower_time = 0.;
	for (long i = 0; i < T; i++) {
		double _t = 0.;
		{
			AutoTimer timer(&_t);
			slower_version(_slots, rnd, context);
		}
		slower_time += _t;
	}
	printf("%.3f + %.3f/per <-> %.3f/per\n", pre_time, faster_time / T, slower_time / T);
	return 0;
}

int main() {
	for (long p : {84961, 82241, 82561, 70913, 84481, 87041})
		run(8192, p);
	return 0;
}

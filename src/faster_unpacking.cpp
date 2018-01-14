#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/NumbTh.h>
#include <HElib/EncryptedArray.h>

#include "CryptGMM/DoublePacking.hpp"
#include "CryptGMM/Timer.hpp"
#include "CryptGMM/HElib.hpp"

#include <iostream>
#include <vector>
void randomize(NTL::ZZX &poly, long p = 3) {
	for (long d = 0; d < NTL::deg(poly); d)
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

int main(int argc, char *argv[]) {
	long m = 8192;
    long p = 70913;
    const long r = 1;
    const long L = 2;
	ArgMapping argmap;
	argmap.arg("m", m, "m");
	argmap.arg("p", p, "p");
	argmap.parse(argc, argv);
    NTL::zz_p::init(p);
    FHEcontext context(m, p, r);
    context.bitsPerLevel = 60;
    buildModChain(context, L);
    std::cerr << "kappa = " << context.securityLevel() << std::endl;
    std::cerr << "slot = " << context.ea->size() << std::endl;
    std::cerr << "degree = " << context.ea->getDegree() << std::endl;

	NTL::ZZX rnd;
	rnd.SetLength(context.zMStar.getPhiM());
	randomize(rnd, p);
	double pre_time = 0.;
	std::vector<GMMPrecompTable> tbls;
	{
		AutoTimer timer(&pre_time);
		tbls = precompute_gmm_tables(context);
	}

	std::vector<long> slots;
	double faster_time = 0.;
	for (long i = 0; i < 10000; i++) {
		double _t;
		{
			AutoTimer timer(&_t);
			faster_version(slots, rnd, tbls, context);
		}
		faster_time += _t;
	}

	std::vector<NTL::zz_pX> _slots;
	double slower_time = 0.;
	for (long i = 0; i < 10000; i++) {
		double _t;
		{
			AutoTimer timer(&_t);
			slower_version(_slots, rnd, context);
		}
		slower_time += _t;
	}
	printf("%.3f + %.3f/per <-> %.3f/per\n", 
		   pre_time, faster_time, slower_time);		
	return 0;
}

#include <HElib/FHE.h>
#include <HElib/FHEContext.h>
#include <HElib/EncryptedArray.h>
#include <HElib/PAlgebra.h>
#include <NTL/ZZX.h>

#include "CryptGMM/HElib.hpp"
long inner_product(NTL::zz_pX const& a, 
                   NTL::zz_pX const& b) 
{
    NTL::zz_p ip(0);
    long deg = NTL::deg(a);
    for (long i = 0; i <= deg; i++) {
        ip += (NTL::coeff(a, i) * NTL::coeff(b, deg - i));
    }
    return ip._zz_p__rep;
}

long inner_product(NTL::ZZX const& a, 
                   NTL::ZZX const& b,
				   long p) 
{
    long deg = NTL::deg(a);
	NTL::ZZ ip(0);
	NTL::ZZ P(p);
    for (long i = 0; i <= deg; i++) {
        ip += NTL::MulMod(NTL::coeff(a, i), NTL::coeff(b, deg - i), P);
    }
    return NTL::to_long(ip) % p;
}

void random_poly(NTL::ZZX &poly, long coeff, long degree)
{
	poly.SetLength(degree);
	for (long i = 0; i < degree; i++)
		NTL::SetCoeff(poly, i, NTL::RandomBnd(coeff) + 1);
}

void test_with_normal_encode() {
    long m = 4096<<2;
    long p = 769;
    FHEcontext context(m, p, 1);
    buildModChain(context, 4);
    FHESecKey sk(context);
    sk.GenSecKey(64);

    const auto &factors = context.alMod.getFactorsOverZZ();
	//auto ea = context.ea;
    EncryptedArray *ea = new EncryptedArray(context, factors[0]);
    const long l = ea->size();
    const long d = ea->getDegree();
	for (long _i = 0; _i < 100; _i++) {
		NTL::ZZX B;
		random_poly(B, p, d);

		std::vector<NTL::ZZX> Vec_A(l);
		std::vector<long> inner_products(l);
		for (long i = 0; i < l; i++) {
			random_poly(Vec_A[i], p, d);
			inner_products[i] = inner_product(Vec_A[i], B, p);
		}

		Ctxt ctx(sk);
		ea->skEncrypt(ctx, sk, std::vector<NTL::ZZX>(l, B));
		NTL::ZZX encoded_A;

		ea->encode(encoded_A, Vec_A);
		ctx.multByConstant(encoded_A);
		std::vector<NTL::ZZX> results;
		ea->decrypt(ctx, sk, results);
		std::vector<long> computed;
		for (auto &s : results) {
			computed.push_back(NTL::to_long(NTL::coeff(s, d - 1)));
		}
		for (long i = 0; i < l; i++) {
			if (computed.at(i) != inner_products.at(i)) {
				std::cout << NTL::deg(Vec_A[i]) << "->";
				std::cout << "computed " << computed.at(i) << " but want " <<
					inner_products.at(i) << std::endl;
			}
		}
	}
}

void test_it() {
    long m = 4096 << 2;
    long p = 769;
    FHEcontext context(m, p, 1);
    buildModChain(context, 4);
    FHESecKey sk(context);
    sk.GenSecKey(64);

	auto ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();

	for (long _i = 0; _i <100; _i++) {
		NTL::zz_p::init(p);

		NTL::zz_pX b;
		NTL::ZZX B;
		NTL::random(b, d);
		NTL::conv(B, b);

		std::vector<NTL::zz_pX> vec_A(l);
		std::vector<long> inner_products(l);
		for (long i = 0; i < l; i++) {
			NTL::random(vec_A[i], d);
			NTL::SetCoeff(vec_A[i], d - 1, 1); // make sure full degree
			inner_products[i] = inner_product(vec_A[i], b);
		}

		NTL::ZZX encoded_A;
		rawEncode(encoded_A, vec_A, context);

		Ctxt ctx(sk);
		sk.Encrypt(ctx, B);
		ctx.multByConstant(encoded_A);

		NTL::ZZX decrypted;
		sk.Decrypt(decrypted, ctx);

		std::vector<NTL::zz_pX> results;
        rawDecode(results, decrypted, context);

		std::vector<long> computed;
		for (auto &s : results) {
			computed.push_back(NTL::coeff(s, d - 1)._zz_p__rep);
		}

		for (long i = 0; i < l; i++) {
			if (computed.at(i) != inner_products.at(i)) {
				std::cout << NTL::deg(vec_A[i]) << "->";
				std::cout << "computed " << computed.at(i) << " but want " <<
					inner_products.at(i) << std::endl;
			}
		}
	}
}

void test_double_pack() {
    long m = 4096;
    long p = 113;
    FHEcontext context(m, p, 1);
    buildModChain(context, 4);
    FHESecKey sk(context);
    sk.GenSecKey(64);

	auto ea = context.ea;
    const long l = ea->size();
    const long d = ea->getDegree();
	const auto &factors = context.alMod.getFactorsOverZZ();

	for (long _i = 0; _i < 1; _i++) {
		NTL::zz_p::init(p);

		std::vector<NTL::zz_pX> vec_A(l);
		std::vector<NTL::zz_pX> vec_B(l);
		for (long i = 0; i < l; i++) {
			NTL::random(vec_B[i], d);
			NTL::SetCoeff(vec_B[i], d - 1, 1); // make sure full degree
        }

		std::vector<long> inner_products(l);
		for (long i = 0; i < l; i++) {
			NTL::random(vec_A[i], d);
			NTL::SetCoeff(vec_A[i], d - 1, 1); // make sure full degree
			long sm = 0;
			for (const auto& vec_b : vec_B) {
				sm += inner_product(vec_A[i], vec_B[i]);
				sm %= p;
			}
			std::cout << sm << " ";
			inner_products[l] = sm;
		}
		std::cout << "\n" << std::endl;

		NTL::ZZX Vec_B;
		Vec_B.SetLength(l * d);
		auto itr = Vec_B.rep.begin();
		for (long i = 0; i < l; i++) {
			long factor = NTL::to_long(factors[i][0]);
			factor = NTL::PowerMod(factor, i, p); // alpha_i^(i+1) mod p
			auto inv = NTL::InvMod(factor, p);
			if (i & 1)
				inv *= -1;
			std::cout << factor << " " << inv << ",";
			for (const auto &b : vec_B[i].rep) {
				NTL::conv(*itr++, inv * b);
			}
		}
		std::cout << '\n';

		NTL::ZZX encoded_A;
		rawEncode(encoded_A, vec_A, context);

		Ctxt ctx(sk);
		sk.Encrypt(ctx, Vec_B);
		ctx.multByConstant(encoded_A);

		NTL::ZZX decrypted;
		sk.Decrypt(decrypted, ctx);

		std::vector<NTL::zz_pX> results;
        rawDecode(results, decrypted, context);

		std::vector<long> computed;
		for (auto &s : results) {
			computed.push_back(NTL::coeff(s, d - 1)._zz_p__rep);
			std::cout << NTL::coeff(s, d - 1) << " ";
		}
		std::cout << "\n";
		// for (long i = 0; i < l; i++) {
		// 	if (computed.at(i) != inner_products.at(i)) {
		// 		std::cout << NTL::deg(vec_A[i]) << "->";
		// 		std::cout << "computed " << computed.at(i) << " but want " <<
		// 			inner_products.at(i) << std::endl;
		// 	}
		// }
	}
}
int main() {
	test_double_pack();
	// auto st = std::clock();
    // test_it();
	// std::cout << (std::clock() - st) / (double)CLOCKS_PER_SEC << std::endl;
    //
	// st = std::clock();
	// test_with_normal_encode();
	// std::cout << (std::clock() - st) / (double)CLOCKS_PER_SEC << std::endl;
	return 0;
}

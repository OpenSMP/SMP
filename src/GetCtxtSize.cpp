#include <iostream>
#include <HElib/FHEContext.h>
#include <HElib/Ctxt.h>
#include <HElib/FHE.h>
void print_size(Ctxt const &ctx) {
	std::stringstream sstream;
	sstream << ctx;
	std::string str = sstream.str();
	std::cout << "|ctxt|" << str.size() / 1024.0 << " KB" << std::endl;
}

int main(int argc, char *argv[]) {
    ArgMapping argmap;
	long m = 8192;
    long p = 8191;
    long L = 2;
	long P = 30;
    argmap.arg("m", m, "m");
    argmap.arg("p", p, "p");
    argmap.arg("P", P, "bits per level");
    argmap.arg("L", L, "L");
    argmap.parse(argc, argv);

	FHEcontext context(m, p, 1);
	context.bitsPerLevel = P;
	buildModChain(context, L);
	FHESecKey sk(context);
	sk.GenSecKey(64);
	Ctxt ctx(sk);
	sk.Encrypt(ctx, NTL::to_ZZX(1));
	print_size(ctx);

	std::cout << "Size after mod-switch to the lowest level: ";
	ctx.modDownToLevel(1);
	print_size(ctx);
	return 0;
}

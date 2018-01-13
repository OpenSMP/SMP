#include "CryptGMM/network/net_io.hpp"
#include <HElib/FHEContext.h>
#include <iostream>
void send_context(std::ostream &s, FHEcontext const& context)
{
    writeContextBase(s, context);
    s << context;
}

FHEcontext receive_context(std::istream &s)
{
	unsigned long m, p, r;
    std::vector<long> gens, ords;
    readContextBase(s, m, p, r, gens, ords);
    FHEcontext context(m, p, r, gens, ords);
    s >> context;
    return context;
}

void receive_context(std::istream &s, FHEcontext **out)
{
	unsigned long m, p, r;
    std::vector<long> gens, ords;
    readContextBase(s, m, p, r, gens, ords);
	*out = new FHEcontext(m, p, r, gens, ords);
    s >> **out;
}

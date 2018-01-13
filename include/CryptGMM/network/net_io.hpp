#ifndef CRYPTGMM_NETWORK_NET_IO_HPP
#define CRYPTGMM_NETWORK_NET_IO_HPP
#include <iosfwd>
class FHEcontext;
void send_context(std::ostream &s, FHEcontext const& context);
FHEcontext receive_context(std::istream &s);
void receive_context(std::istream &s, FHEcontext **out);
#endif // CRYPTGMM_NETWORK_NET_IO_HPP

#ifndef CRYPTGMM_NETWORK_NET_IO_HPP
#define CRYPTGMM_NETWORK_NET_IO_HPP
class FHEContext;
class Ctxt;
template <class Imp>
class NetIO {
public:
    ~NetIO() {}

    void close() {
        derived().close();
    }

    virtual void send_context(FHEContext const& context) = 0;

    virtual void receive_context(FHEContext &context) = 0;
    
    virtual void send_ctxt(Ctxt const& ctxt) = 0;

    virtual void receive_ctxt(Ctxt &ctxt) = 0;
private:
    Imp& derived() {
        return *static_cast<Imp*>(this);
    }
};
#endif // CRYPTGMM_NETWORK_NET_IO_HPP

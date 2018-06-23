#include <NTL/lzz_pX.h>
#include <HElib/FHE.h>
#include <HElib/EncryptedArray.h>
#include <vector>
#include <memory>
using RX = NTL::zz_pX;
using ModRX = NTL::zz_pXModulus;

class BasicAutomorphPrecon {
  Ctxt ctxt;
  NTL::xdouble noise;
  std::vector<DoubleCRT> polyDigits;

public:
  BasicAutomorphPrecon(const Ctxt& _ctxt) : ctxt(_ctxt), noise(1.0)
  {
    FHE_TIMER_START;
    if (ctxt.parts.size() >= 1) assert(ctxt.parts[0].skHandle.isOne());
    if (ctxt.parts.size() <= 1) return; // nothing to do

    ctxt.cleanUp();
    const FHEcontext& context = ctxt.getContext();
    const FHEPubKey& pubKey = ctxt.getPubKey();
    long keyID = ctxt.getKeyID();

    // The call to cleanUp() should ensure that this assertions passes.
    assert(ctxt.inCanonicalForm(keyID));

    // Compute the number of digits that we need and the esitmated
    // added noise from switching this ciphertext.
    long nDigits;
    std::tie(nDigits, noise)
      = ctxt.computeKSNoise(1, pubKey.keySWlist().at(0).ptxtSpace);

    double logProd = context.logOfProduct(context.specialPrimes);
    noise += ctxt.getNoiseVar() * xexp(2*logProd);

    // Break the ciphertext part into digits, if needed, and scale up these
    // digits using the special primes.

    ctxt.parts[1].breakIntoDigits(polyDigits, nDigits);
  }

  
  shared_ptr<Ctxt>
  automorph(long k) const
  {
    FHE_TIMER_START;

    // A hack: record this automorphism rather than actually performing it
    if (isSetAutomorphVals()) { // defined in NumbTh.h
      recordAutomorphVal(k);
      return make_shared<Ctxt>(ctxt);
    }

    if (k==1 || ctxt.isEmpty()) return make_shared<Ctxt>(ctxt);// nothing to do

    const FHEcontext& context = ctxt.getContext();
    const FHEPubKey& pubKey = ctxt.getPubKey();
    shared_ptr<Ctxt> result = make_shared<Ctxt>(ZeroCtxtLike, ctxt); // empty ctxt
    result->noiseVar = noise; // noise estimate

    if (ctxt.parts.size()==1) { // only constant part, no need to key-switch
      CtxtPart tmpPart = ctxt.parts[0];
      tmpPart.automorph(k);
      tmpPart.addPrimesAndScale(context.specialPrimes);
      result->addPart(tmpPart, /*matchPrimeSet=*/true);
      return result;
    }

    // Ensure that we have a key-switching matrices for this automorphism
    long keyID = ctxt.getKeyID();
    if (!pubKey.isReachable(k,keyID)) {
      throw std::logic_error("no key-switching matrices for k="+std::to_string(k)
                             + ", keyID="+std::to_string(keyID));
    }

    // Get the first key-switching matrix for this automorphism
    const KeySwitch& W = pubKey.getNextKSWmatrix(k,keyID);
    long amt = W.fromKey.getPowerOfX();

    // Start by rotating the constant part, no need to key-switch it
    CtxtPart tmpPart = ctxt.parts[0];
    tmpPart.automorph(amt);
    tmpPart.addPrimesAndScale(context.specialPrimes);
    result->addPart(tmpPart, /*matchPrimeSet=*/true);

    // Then rotate the digits and key-switch them
    vector<DoubleCRT> tmpDigits = polyDigits;
    for (auto&& tmp: tmpDigits) // rotate each of the digits
      tmp.automorph(amt);

    result->keySwitchDigits(W, tmpDigits); // key-switch the digits

    long m = context.zMStar.getM();
    if ((amt-k)%m != 0) { // amt != k (mod m), more automorphisms to do
      k = MulMod(k, InvMod(amt,m), m); // k *= amt^{-1} mod m
      result->smartAutomorph(k);       // call usual smartAutomorph
    }
    return result;
  }
};

void frob(RX& f, long j, long t, ModRX const& mod)
{
    long d = NTL::deg(mod);
    j = mcMod(j, d);
    RX H = PowerMod(RX(1, 1), NTL::power_ZZ(t, j), mod);
    f = NTL::CompMod(f, H, mod);
}

struct ExtractAux {
    long m;
    long t;
    long d;
};

void extract(Ctxt &ctx, 
             RX const& alpha,
             RX const& u,
             RX const& action,
             ModRX const& G,
             const EncryptedArray *ea,
             ExtractAux const& aux)
{
    const long l = ea->size();
    RX alpha_u = NTL::MulMod(alpha, u, G);
    Ctxt ans(ZeroCtxtLike, ctx);
    BasicAutomorphPrecon precon(ctx);
    RX C;
    zzX C_p;
    for (long j = 0; j < aux.d; ++j) {
        long amt = NTL::PowerMod(aux.t, mcMod(j, aux.d), aux.m);
        auto frb_ctx = precon.automorph(amt);

        auto power_j(alpha_u);
        if (j > 0)
            frob(power_j, j, aux.t, G);

        NTL::MulMod(C, power_j, action, G);
        convert(C_p, C);

        std::vector<zzX> slots(l, C_p);
        ea->encode(C_p, slots); // This part took most of the computation time
        frb_ctx->multByConstant(C_p);
        ans += *frb_ctx;
    }
    std::swap(ctx, ans);
}

int main() {
    const long m = 8192;
    const long t = 67073;
    NTL::zz_p::init(t);
    FHEcontext context(m, t, 1);
    buildModChain(context, 3);
    FHESecKey sk(context);
    sk.GenSecKey(64);
    addFrbMatrices(sk);
    FHEPubKey const& pk{sk};
    const EncryptedArray *ea = context.ea;
    const long d = ea->getDegree();
    const long l = ea->size();
    std::cout << "d = " << d << "\n";

    const auto &factors = context.alMod.getFactorsOverZZ();
    std::vector<ModRX> modFt(l);
    ModRX modPhiX;
    for (long i = 0; i < l; ++i) {
        NTL::build(modFt[i], NTL::conv<RX>(factors[i]));
    }
    NTL::build(modPhiX, NTL::conv<RX>(context.zMStar.getPhimX()));

    NTL::ZZX msg;
    RX msg_p;
    {
        NTL::ZZX _msg;
        _msg.SetLength(d);
        for (long i = 0; i < d; ++i) _msg[i] = i + 1;
        std::vector<NTL::ZZX> slots(l, _msg);
        NTL::conv(msg_p, _msg);
        ea->encode(msg, slots);
    }
    Ctxt ctx(sk);
    sk.Encrypt(ctx, msg);
    NTL::ZZX alpha;
    RX alpha_p;
    {
        std::vector<NTL::ZZX> Ls(d, NTL::ZZX(0));
        Ls.back() = NTL::ZZX(1);
        std::vector<NTL::ZZX> coeffs;
        ea->buildLinPolyCoeffs(coeffs, Ls);
        std::vector<NTL::ZZX> slots(l, coeffs.front());
        ea->encode(alpha, slots);
        NTL::conv(alpha_p, coeffs.front());
    }

    RX u = msg_p;
    std::vector<Ctxt> ctxs(d, ctx);
    auto _st = std::clock();
    for (long i = 0; i < d; ++i) {
        RX action_p = NTL::conv<RX>(NTL::ZZX(i, 1));
        extract(ctxs[i], alpha_p, u, action_p, modFt[0], ea, {m, t, d});
        if (i > 0)
            ctxs[0] += ctxs[i];
    }
    auto _end = std::clock();

    {
        std::vector<NTL::ZZX> slots;
        ea->decrypt(ctxs[0], sk, slots);
        std::cout << slots[0] << std::endl;
    }

    std::cout << "kappa " << context.securityLevel() << std::endl;
    std::cout << ctxs[0].log_of_ratio() << std::endl;
    std::cout << "extract " << (_end - _st) / (double) CLOCKS_PER_SEC << std::endl;
    printAllTimers(std::cout);
    return 0;
}

//
// Created by riku on 6/20/18.
//
#include "HElib/FHE.h"
#include "HElib/FHEContext.h"
#include "HElib/PAlgebra.h"
#include "HElib/EncryptedArray.h"
#include <NTL/ZZX.h>

#include "SMP/HElib.hpp"

/*
void mapToSlots() {
    for (long i = 0; i < nSlots; i++) {
        long t = zMStar.ith_rep(i);
        long tInv = InvMod(t, m);
        RX ct_rep;
        PowerXMod(ct_rep, tInv, G);
        RE ct;
        conv(ct, ct_rep);
        REX Qi;
        SetCoeff(Qi, 1, 1);
        SetCoeff(Qi, 0, -ct);
        mappingData.rmaps[i] = Qi;
    }
}
*/

int main() {
    FHEcontext context(8192, 67073, 1);

    PA_zz_p::R::init(context.zMStar.getP());
    MappingData<PA_zz_p> mappingData;
    PA_zz_p::RX G;
    NTL::conv(G, context.alMod.getFactorsOverZZ()[0]);
    context.alMod.getDerived(PA_zz_p()).mapToSlots(mappingData, G);

    const EncryptedArray *ea = context.ea;
    std::vector<PA_zz_p::RX> slots(ea->size());
    for (int i = 0; i < ea->size(); ++i) {
        slots[i].SetLength(ea->getDegree());
        for (int j = 0; j < ea->getDegree(); j++) {
            NTL::SetCoeff(slots[i], j, NTL::RandomBnd(67073));
        }
    }

    NTL::ZZX encoded;
    std::vector<NTL::ZZX> _slots;
    convert(_slots, slots);

    ea->encode(encoded, _slots);

    NTL::ZZX _encoded;
    rawEncode(_encoded, slots, context);

    if (encoded == _encoded)
        std::cout << "eq\n";
    else
        std::cout << "ne\n";
    return 0;
}


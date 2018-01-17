//
// Created by Lu WJ on 5/7/2017 AD.
//

#ifndef CRYPTGMM_DOUBLEPACKING_HPP
#define CRYPTGMM_DOUBLEPACKING_HPP
#include "SMP/Matrix.hpp"
#include "NTL/ZZX.h"
#include <vector>

class FHESecKey;
class EncryptedArray;
namespace internal {
    struct BlockId {
        int x, y; 
    };

    struct PackedRows {
        std::vector<NTL::zz_pX> polys;
        long num_duplication;
    };

    PackedRows partition(Matrix const& matrix, 
                         BlockId const& blk, 
						 const EncryptedArray *packer,
                         const bool backward);
}
#endif //CRYPTGMM_DOUBLEPACKING_HPP

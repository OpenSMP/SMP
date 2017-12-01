//
// Created by Lu WJ on 5/7/2017 AD.
//

#ifndef CRYPTGMM_DOUBLEPACKING_HPP
#define CRYPTGMM_DOUBLEPACKING_HPP
#include "CryptGMM/Matrix.hpp"
#include "NTL/ZZX.h"
#include <vector>

class FHESecKey;
using Packer = EncryptedArray;
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
                         Packer const& packer,
                         const bool backward);
}
#include "CryptGMM/internal/DoublePacking.hpp"
#endif //CRYPTGMM_DOUBLEPACKING_HPP

//
// Created by Lu WJ on 5/7/2017 AD.
//

#ifndef CRYPTGMM_DOUBLEPACKING_HPP
#define CRYPTGMM_DOUBLEPACKING_HPP
#include "CryptGMM/CowPtr.hpp"
#include "CryptGMM/Matrix.hpp"
#include "NTL/ZZX.h"
#include <vector>

class Ctxt;
class FHEPubKey;
class FHESecKey;
class EncryptedArray;
class GMMController;

typedef EncryptedArray Packer;
typedef CowPtr<Ctxt> Cipher;
namespace internal {
    struct BlockId {
        int x, y; 
    };

    struct PackedRows {
        std::vector<NTL::ZZX> polys;
        std::vector<size_t> row_indices;
    };

    PackedRows partition(Matrix const& matrix, 
                         BlockId const& blk, 
                         Packer const& packer,
                         const bool backward);
}

class DoublePackedMat {
public:
    DoublePackedMat();

    ~DoublePackedMat();

    void pack(Matrix const& matrix, Packer const& packer);

    void encrypt(Matrix const&matrix, FHEPubKey const& key);

    void encrypt(Matrix const&matrix, FHESecKey const& key);

    long num_rows() const;

    long num_cols() const;

    bool load(std::istream &stream);

    bool save(std::ostream &stream) const;

    friend class GMMController;
protected:
    class Imp;
    CowPtr<Imp> imp_;
};

#endif //CRYPTGMM_DOUBLEPACKING_HPP

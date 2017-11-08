#ifndef CRYPTGMM_GMMCONTROLLER_HPP
#define CRYPTGMM_GMMCONTROLLER_HPP
#include "CowPtr.hpp"
#include "Matrix.hpp"
#include <vector>

class Ctxt;
class FHESecKey;
class FHEPubKey;
typedef CowPtr<Ctxt> Cipher;

struct GMMResult {
    enum class Status {
        OK, 
        BAD_ARG,
    };
    std::vector<Cipher> ciphers;
    long num_cols;
    long num_rows;
    Status status;
};

class DoublePackedMat;
class GMMController {
public:
    GMMController() {}

    ~GMMController() {}
    /// @return mat1 * mat2
    GMMResult multiply(DoublePackedMat const& mat1,
                       DoublePackedMat const& mat2) const;

    void decrypt(PlainVec *out, GMMResult const& result, FHESecKey const& key) const;
};

#endif

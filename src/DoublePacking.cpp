#include "CryptGMM/DoublePacking.hpp"
#include "CryptGMM/GMMController.hpp"
#include "CryptGMM/internal/DoublePacking.hpp"

#include "HElib/FHE.h"
#include "HElib/Ctxt.h"

#include <list>
class DoublePackedMat::Imp {
public:
    typedef struct {
        std::vector<Cipher> parts;
    } CipherBlock;

    Imp() : num_rows_(0), num_cols_(0) {}

    ~Imp() {}

    void pack(Matrix const& matrix, Packer const& packer, bool backward = false) {
        num_rows_ = matrix.NumRows();
        num_cols_ = matrix.NumCols();
        packed_.clear();

        const long d = packer.getDegree();
        const long l = packer.size();
        const long MAX_X = (num_rows_ + d - 1) / d;
        const long MAX_Y = (num_cols_ + l - 1) / l;
        for (int y = 0; y < MAX_Y; y++) {
            std::vector<NTL::ZZX> polys;
            for (int x = 0; x < MAX_X; x++) {
                internal::BlockId blk{x, y};
                std::vector<NTL::ZZX> slots = internal::partition(matrix, blk, packer, backward);
                NTL::ZZX poly;
                packer.encode(poly, slots);
                polys.emplace_back(poly);
            }
            packed_.emplace_back(polys);
        }
    }

    template <class Key>
    void encrypt(Matrix const& matrix, Key const& key) {
        const EncryptedArray *packer = key.getContext().ea;
        pack(matrix, *packer, /*backward =*/true);
        encrypted_.clear();
        for (const auto& packed : packed_) {
            CipherBlock block;
            for (const auto& poly : packed) {
                Ctxt *ctxt = new Ctxt(key);
                key.Encrypt(*ctxt, poly);
                block.parts.emplace_back(ctxt);
            }
        }
        /// clear up plain values
        packed_.clear();
    }

    bool save(std::ostream &stream) const { return false; }

    bool load(std::istream &stream) { return true; }

    long num_cols() const { return num_cols_; }

    long num_rows() const { return num_rows_; }

    const std::vector<CipherBlock>& encrypted_content() const { return encrypted_; }

    bool is_encrypted() const { return !encrypted_.empty(); }

    bool is_packed() const { return !packed_.empty(); }

private:
    long num_rows_;
    long num_cols_;
    std::vector<CipherBlock> encrypted_;
    std::list<std::vector<NTL::ZZX>> packed_;
};

// GMMResult 
// GMMController::multiply(DoublePackedMat const& mat1,
//                         DoublePackedMat const& mat2) const {
//     GMMResult result;
//     if (mat1.type_ != DoublePackedMat::Type::LEFT &&
//         mat2.type_ != DoublePackedMat::Type::RIGHT) {
//         result.status = GMMResult::Status::BAD_ARG;
//         return result;
//     }
//     const auto& enc_mat1 = mat1.imp_->encrypted_content();
//     const auto& enc_mat2 = mat2.imp_->encrypted_content();
//     
//     size_t num_rows = enc_mat1.size();
//     size_t num_blks = enc_mat2.size();
//     for (size_t r = 0; r < num_rows; r++) {
//         auto cpy(enc_mat1[r]);
//         for (size_t b = 0; b < num_blks; b++) {
//             assert(cpy.parts.size() == mat2[b].parts.size());
//         }
//     }
// }
//
void DoublePackedMat::pack(Matrix const& matrix, Packer const& packer) {
    imp_->pack(matrix, packer);
}

void DoublePackedMat::encrypt(Matrix const& mat, FHEPubKey const& key) {
    imp_->encrypt(mat, key);
}

void DoublePackedMat::encrypt(Matrix const& mat, FHESecKey const& key) {
    imp_->encrypt(mat, key);
}

long DoublePackedMat::num_rows() const {
    return imp_->num_rows();
}

long DoublePackedMat::num_cols() const {
    return imp_->num_cols();
}

bool DoublePackedMat::load(std::istream &stream) {
    return imp_->load(stream);
}

bool DoublePackedMat::save(std::ostream &stream) const {
    return imp_->save(stream);
}

DoublePackedMat::DoublePackedMat() {
    imp_ = CowPtr<Imp>(new Imp());
}

DoublePackedMat::~DoublePackedMat() { }

//
// Created by Lu WJ on 5/6/2017 AD.
//

#ifndef CRYPTCONV_ENCVEC_HPP
#define CRYPTCONV_ENCVEC_HPP
#include "Tensor.hpp"
#include "CowPtr.hpp"
#include <iosfwd>

class FHEPubKey;
class FHESecKey;
class EncVec {
public:
    EncVec();

    ~EncVec();

    void encrypt(const PlainVec &vec, const FHEPubKey &key);

    void decrypt(PlainVec *out, const FHESecKey &key) const;

    EncVec &operator+=(const EncVec &oth);

    EncVec &mul(const long scalar);

    EncVec &rotate(const long offset);

    size_t num_ctxts() const;

    long length() const;

    long cyclotomic() const;

    bool save_to(const std::string &file) const;

    bool save_to(std::ostream &os) const;

    bool load_from(std::istream &is, const FHEPubKey &key);

    bool load_from(const std::string &file, const FHEPubKey &key);
private:
    class Imp;
    typedef CowPtr<Imp> ImpPtr;
    ImpPtr imp_;
};
#endif //CRYPTCONV_ENCVEC_HPP

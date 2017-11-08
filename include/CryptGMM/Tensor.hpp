//
// Created by Lu WJ on 4/25/2017 AD.
//

#ifndef CRYPTCONV_TENSOR_HPP
#define CRYPTCONV_TENSOR_HPP
#include <vector>
#include <iosfwd>
#include "Matrix.hpp"
class FHEPubKey;
class FHESecKey;
long conv_size(long n, long h, long s);

class Tensor {
public:
    Tensor(const std::string &file);

    Tensor(long h, long w, long c);

    Tensor(const std::vector<Matrix> &init);

    Tensor(const Tensor &oth) : h_(oth.h_), w_(oth.w_), c_(oth.c_),
                                raw_(oth.raw_) {}

    ~Tensor();

    void relu();

    bool operator==(const Tensor &oth) const;

    bool operator!=(const Tensor &oth) const { return ! (*this == oth); }

    void set(long i, long j, long k, const val_t &val);

    val_t get(long i, long j, long k) const;

    Matrix conv2d(const Tensor &kernel, const long stride) const;

    Tensor conv(const std::vector<Tensor> &kernels, const long stride) const;

    Tensor depthwise_conv(const Tensor &kernels, const long stride) const;

    Tensor pointwise_conv(const std::vector<Tensor> &kernels) const;

    long channels() const;

    long rows() const;

    long columns() const;

    bool save(const std::string &save_to) const;

    friend std::ostream& operator <<(std::ostream &out, const Tensor &);

    const Matrix& channel(long c) const;

    Matrix conv2d(const Matrix &channel, const Matrix &kernel, const long stride) const;

    void reshape(const PlainVec &vec);
private:
    Matrix pointwise_conv(const Tensor &kernel) const;

private:
    bool parse_header(const std::string &file);

private:
    long h_, w_, c_;
    std::vector<Matrix> raw_;
};

bool is_same(const Matrix &a, const Matrix &b);
#endif //CRYPTCONV_TENSOR_HPP

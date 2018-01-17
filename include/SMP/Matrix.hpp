//
// Created by Lu WJ on 5/7/2017 AD.
//

#ifndef CRYPTCONV_MATRIX_HPP
#define CRYPTCONV_MATRIX_HPP
#include <NTL/matrix.h>
typedef long val_t;
typedef NTL::Mat<val_t> Matrix;
typedef NTL::Vec<val_t> PlainVec;
void randomize(Matrix *mat);
void zeros(Matrix *mat, long N);
void add(Matrix *mat, const Matrix &b);
PlainVec mul(const Matrix &a, const PlainVec &b);
Matrix mul(const Matrix &a, const Matrix &b);
void add(PlainVec *a, const PlainVec &b);
void transpose(Matrix *T, const Matrix &mat);
bool is_same(const Matrix &a, const Matrix &b);
bool is_same(const Matrix &a, const Matrix &b, long modulus);
bool load_matrix(Matrix *mat, const std::string &file);
bool load_matrix(std::istream& in, Matrix *mat);
bool save_matrix(const Matrix &mat, const std::string &file);
bool save_matrix(std::ostream &out, const Matrix &mat);
bool load_vector(PlainVec *vec, const std::string &file);
bool load_vector(std::istream &in, PlainVec *vec);
bool save_vector(const PlainVec &vec, const std::string &file);
bool save_vector(std::ostream &out, const PlainVec &vec);
#endif //CRYPTCONV_MATRIX_HPP

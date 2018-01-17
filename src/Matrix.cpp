//
// Created by Lu WJ on 5/7/2017 AD.
//
#include "SMP/Matrix.hpp"
#include "SMP/literal.hpp"

#include <NTL/ZZ.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <fstream>

void zeros(Matrix *mat, long N) {
    if (!mat) return;
    mat->SetDims(N, N);
    for (long i = 0; i < N; i++) {
        for (long j = 0; j < N; j++)
            (*mat)[i][j] = 0L;
    }
}

void add(PlainVec *a, const PlainVec &b) {
    if (!a or a->length() != b.length())
        return;
    for (long i = 0; i < b.length(); i++)
        a->put(i, a->get(i) + b[i]);
}

void add(Matrix *a, const Matrix &b) {
    if (!a) return;
    assert(a->NumRows() == b.NumRows());
    assert(a->NumCols() == b.NumCols());
    for (long r = 0; r < b.NumRows(); r++) {
        for (long c = 0; c < b.NumCols(); c++) {
            (*a)[r][c] += b[r][c];
        }
    }
}

void transpose(Matrix *T, const Matrix &mat) {
    if (!T) return;
    T->kill();
    T->SetDims(mat.NumCols(), mat.NumRows());
    for (long r = 0; r < mat.NumRows(); r++) {
        for (long c = 0; c < mat.NumCols(); c++)
            (*T)[c][r] = mat[r][c];
    }
}

void randomize(Matrix *mat) {
    if (!mat)
        return;
    for (long r = 0; r < mat->NumRows(); r++) {
        for (long c = 0; c < mat->NumCols(); c++) {
            (*mat)[r][c] = NTL::RandomBnd(20) - 10;
        }
    }
}

bool is_same(const Matrix &a, const Matrix &b) {
    if (a.NumRows() != b.NumRows() || a.NumCols() != b.NumCols())
        return false;
    for (long r = 0; r < a.NumRows(); r++) {
        for (long c = 0; c < a.NumCols(); c++) {
            if (a[r][c] != b[r][c]) {
                std::cerr << a[r][c] << "!=" << b[r][c] << "(" << r << "," << c << ")" << std::endl;
                return false;
            }
        }
    }
    return true;
}

bool is_same(const Matrix &a, const Matrix &b, long modulus) {
    if (a.NumRows() != b.NumRows() || a.NumCols() != b.NumCols())
        return false;
    for (long r = 0; r < a.NumRows(); r++) {
        for (long c = 0; c < a.NumCols(); c++) {
            long diff = (a[r][c] - b[r][c]) % modulus;
            if (diff != 0) {
                std::cerr << a[r][c] << "!=" << b[r][c] << "(" << r << "," << c << ")" << std::endl;
                return false;
            }
        }
    }
    return true;
}

PlainVec mul(const Matrix &m, const PlainVec &v) {
    assert(m.NumCols() == v.length());
    PlainVec result;
    result.SetLength(m.NumRows());
    for (long r = 0; r < m.NumRows(); r++) {
        val_t sum{0};
        for (long c = 0; c < m.NumCols(); c++) {
            sum += m[r][c] * v[c];
        }
        result.put(r, sum);
    }
    return result;
}

Matrix mul(const Matrix &a, const Matrix &b) {
    assert(a.NumCols() == b.NumRows());
    Matrix result;
    result.SetDims(a.NumRows(), b.NumCols());
    for (long r = 0; r < a.NumRows(); r++) {
        for (long c = 0; c < b.NumCols(); c++) {
            val_t sum{0};
            for (long k = 0; k < b.NumRows(); k++)
                sum += a[r][k] * b[k][c];
            result[r][c] = sum;
        }
    }
    return result;
}

bool load_matrix(std::istream &in, Matrix *mat) {
    if (in.bad() || in.eof() || !mat)
        return false;

    std::string line;
    int line_num = 0;
    for (long r = 0; !in.eof() && r < mat->NumRows(); r++) {
        std::getline(in, line);
        auto fields = splitBySpace(line);
        if (fields.size() != mat->NumCols()) {
            std::cerr << "Warn: need " << mat->NumCols() << " columns, but got " << fields.size() << 
                " values in the " << line_num + 1<< " line." << std::endl;
            return false;
        }
        size_t pos;
        for (size_t c = 0; c < fields.size(); c++) {
            long val = std::stol(fields[c], &pos, 10);
            if (pos != fields[c].size()) {
                std::cerr << "Warn: invalid value " << fields[c] << " in line " << line_num + 1 << std::endl;
                val = 0;
            }
            (*mat)[line_num][c] = val;
        }
        line_num += 1;
    }

    if (line_num != mat->NumRows()) {
        std::cerr << "Warn: need " << mat->NumRows() << " rows but got " << line_num << " rows." << std::endl;
        return false;
    }
    return true;
}

bool load_matrix(Matrix *mat, const std::string &file) {
    if (!mat) {
        std::cerr << "Error: can not load matrix in an empty pointer." << std::endl;
        return false;
    }

    std::ifstream in(file);
    if (!in.is_open()) {
        std::cerr << "Error: can not open file " << file << std::endl;
        return false;
    }

    std::string header_line;
    std::getline(in, header_line);
    std::stringstream header(header_line);
    char sharp;
    long rows, cols;
    header >> sharp >> rows >> cols;
    if (sharp != '#' || rows <= 0 || cols <= 0) {
        std::cerr << "Error: invalid header " << file << std::endl;
        return false;
    }
    mat->SetDims(rows, cols);
    bool ok = load_matrix(in, mat);
    in.close();
    return ok;
}

bool save_matrix(std::ostream &out, const Matrix &mat) {
    if (out.bad() || out.eof())
        return false;
    if (mat.NumRows() <= 0 || mat.NumCols() <= 0)
        return false;
    for (long r = 0; r < mat.NumRows(); r++) {
        for (long c = 0; c + 1 < mat.NumCols(); c++) {
            out << mat[r][c] << " ";
        }
        out << mat[r][mat.NumCols() - 1] << std::endl;
    }
    return out.good();
}

bool save_matrix(const Matrix &mat, const std::string &file) {
    std::ofstream out(file);
    if (!out.is_open()) {
        std::cerr << "Error: can not open file " << file << std::endl;
        return false;
    }
    out << "#" << mat.NumRows() << " " << mat.NumCols() << std::endl; 
    bool ok = save_matrix(out, mat);
    out.close();
    return ok;
}

bool load_vector(std::istream &in, PlainVec *vec) {
    if (in.bad() || in.eof() || !vec)
        return false;
    std::string line;
    std::getline(in, line);
    const std::vector<std::string> fields = splitBySpace(line);
    if (fields.size() != vec->length()) {
        std::cerr << "Error: need " << vec->length()
                  << " elements but get " << fields.size() << "." << std::endl;
        return false;
    }
    size_t pos;
    for (size_t c = 0; c < fields.size(); c++) {
        long val = std::stol(fields.at(c), &pos, 10);
        if (pos != fields[c].size()) {
            std::cerr << "Warn: invalid value: " << fields[c] << std::endl;
            val = 0L;
        }
        vec->put(c, val);
    }
    return true;
}

bool load_vector(PlainVec *vec, const std::string &file) {
    std::ifstream in(file);
    if (!in.is_open()) {
        std::cerr << "Error: can not open file " << file << std::endl;
        return false;
    }
    std::string header;
    std::getline(in, header);
    std::stringstream sstream(header);
    char sharp = '\0';
    long length = 0;
    sstream >> sharp >> length;
    if (sharp != '#' || length <= 0) {
        std::cerr << "Error: invalid header of " << file << std::endl;
        return false;
    }

    vec->SetLength(length);
    bool ok = load_vector(in, vec);
    in.close();
    return ok;
}

bool save_vector(const PlainVec &vec, const std::string &file) {
    std::ofstream out(file);
    if (!out.is_open()) {
        std::cerr << "Error: can not open file " << file << std::endl;
        return false;
    }

    out << "#" << vec.length() << std::endl;
    bool ok = save_vector(out, vec);
    out.close();
    return ok;
}

bool save_vector(std::ostream &out, const PlainVec &vec) {
    if (out.bad())
        return false;
    if (vec.length() == 0)
        return true;
    for (long c = 0; c + 1< vec.length(); c++)
        out << vec[c] << " ";
    out << vec[vec.length() - 1] << std::endl;
    return out.good();
}

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once

#include "sys.h"
#include <map>

namespace nla3d {
namespace math {
// NOTE: SparceMatrix implementation has indexes started from 1 !

// TODO:
// SparseMatrix in nla3d is primarly used to provide sparse matrix data
// to MKL's PARDISO solver. That's the reason to use uint32 type for indexing.
// It's more general to use size_t for indexes. To manage with it there are several ways:
// 1. Use templated calsses
// 2. Use typedef int SparseIndex_t
// 3. Use #define SparseIndex_t int


class SparsityInfo {
  public:
    SparsityInfo(uint32 _nrows, uint32 _ncols, uint32 _max_in_row);

    // add information than (_row, _column) element has non-zero value
    // row and column positions are started from 1
    void add(uint32 _i, uint32 _j);
    // perform compression procedure. After that positions of non-zero elements can't be changed
    void compress();
    bool isCompressed();

    uint32 nElementsInRow(uint32 _row);
    uint32 getIndex(uint32 _i, uint32 _j);

    friend class BaseSparseMatrix;
    friend class SparseMatrix;
    friend class SparseSymMatrix;
  private:

    // Data arrays to implement compressed sparse row format with 3 arrays (3-array CSR).
    // Format was implemented by using MKL manual.
    // Here we don't need to know excactly values. We need to store only info about non-zeros
    // elemenets in matrix.
    // columns arrays stores column number of element: columns[element_number] = element_column
    uint32* columns = nullptr;
    // index of first element in row: iofeir[row-1] first element in row,
    // ifeir[row]-1 - last element in row
    uint32* iofeir = nullptr;
    // number of non-zero elements in the matrix. Also, size of values array.
    // If numberOfValues > 0 that means that the matrix is ready for work (a training has been completed already)
    uint32 numberOfValues = 0;
    uint32 maxInRow = 0;
    bool compressed = false;

    // size of the matrix
    uint32 nRows;
    uint32 nColumns;

    static const uint32 invalid;
};



// general Sparse Matrix in a row-oriented data format
// The implementation of compressed sparse row format with 3 arrays (3-array CSR).
class BaseSparseMatrix {
  public:
    BaseSparseMatrix(uint32 nrows, uint32 ncolumns, uint32 max_in_row = 100);
    BaseSparseMatrix(std::shared_ptr<SparsityInfo> spar_info);
    ~BaseSparseMatrix();
    
    void compress();

    void add(uint32 _i, uint32 _j);
    void addValue (uint32 _i, uint32 _j, double value);

    void clear();

    // debug methods
    void printInternalData (std::ostream& out);

    void zero();

    // getters
    double* getValuesArray();
    uint32* getColumnsArray();
    uint32* getIofeirArray();
    uint32 getNumberOfValues();
    uint32 getNumberOfRows();
    uint32 getNumberOfColumns();
    std::shared_ptr<SparsityInfo> getSparsityInfo();

    bool isCompressed();

  protected:
    // Data arrays to implement compressed sparse row format with 3 arrays (3-array CSR).
    double *values = nullptr;
    // is our matrix is symmetric. If so only upper triangle is stored

    std::shared_ptr<SparsityInfo> si;
};


class SparseMatrix : public BaseSparseMatrix {
  public:
    SparseMatrix(uint32 nrows, uint32 ncolumns, uint32 max_in_row = 100);
    SparseMatrix(std::shared_ptr<SparsityInfo> spar_info);
    
    void add(uint32 _i, uint32 _j);
    void addValue (uint32 _i, uint32 _j, double value);


    // TODO: It would be better to use a BLAS routine for sparse matrices for speedup
    double mult_vec_i (double *vec, uint32 i);
    double transpose_mult_vec_i (double *vec, uint32 i);
    void transpose_mult_vec (double *vec, double *res);

    // debug methods
    void print (std::ostream& out);

    // getters
    // return reference to the value, if the value doesn't exist in the matrix - fatal error
    double& operator()(uint32 _i, uint32 _j);

    // return the value, if the value doesn't exists - return 0.0
    double value(uint32 _i, uint32 _j) const;
};


class SparseSymMatrix : public BaseSparseMatrix {
  public:
    SparseSymMatrix(uint32 nrows, uint32 max_in_row = 100);
    SparseSymMatrix(std::shared_ptr<SparsityInfo> spar_info);
    
    void add(uint32 _i, uint32 _j);
    void addValue (uint32 _i, uint32 _j, double value);


    double mult_vec_i(double *vec, uint32 i);

    // debug methods
    void print (std::ostream& out);

    // getters
    // return reverence to the value, if the value doesn't exist in the matrix - fatal error
    double& operator()(uint32 _i, uint32 _j);

    // return the value, if the value doesn't exists - return 0.0
    double value(uint32 _i, uint32 _j) const;
};


inline bool SparsityInfo::isCompressed() {
    return compressed;
}

inline uint32 SparsityInfo::nElementsInRow(uint32 _row) {
    assert(iofeir != nullptr);
    assert(_row > 0 && _row <= nRows);
    return iofeir[_row] - iofeir[_row-1];
}


inline void BaseSparseMatrix::zero() {
  assert(si);
  assert(values);
  std::fill_n(&values[0], si->numberOfValues, 0.0);
}

inline double* BaseSparseMatrix::getValuesArray() {
  assert(values != nullptr);
  return values;
}

inline uint32* BaseSparseMatrix::getColumnsArray() {
  assert(si);
  assert(si->compressed);
  return si->columns;
}

inline uint32* BaseSparseMatrix::getIofeirArray() {
  assert(si);
  assert(si->compressed);
  return si->iofeir;
}

inline uint32 BaseSparseMatrix::getNumberOfValues() {
  assert(si);
  return si->numberOfValues;
}

inline uint32 BaseSparseMatrix::getNumberOfRows() {
  assert(si);
  return si->nRows;
}

inline uint32 BaseSparseMatrix::getNumberOfColumns() {
  assert(si);
  return si->nColumns;
}

inline std::shared_ptr<SparsityInfo> BaseSparseMatrix::getSparsityInfo() {
  assert(si);
  return si;
}

inline bool BaseSparseMatrix::isCompressed() {
  assert(si);
  return si->compressed;
}


inline void SparseMatrix::add(uint32 _i, uint32 _j) {
  assert(si);
  si->add(_i, _j);
}

inline void SparseMatrix::addValue(uint32 _i, uint32 _j, double value) {
  this->operator()(_i, _j) += value;
}


inline void SparseSymMatrix::add(uint32 _i, uint32 _j) {
  assert(si);

  // ensure that we work in upper triangle
	if (_i > _j) std::swap(_i, _j);

  si->add(_i, _j);
}

inline void SparseSymMatrix::addValue(uint32 _i, uint32 _j, double value) {
  this->operator()(_i, _j) += value;
}

 
} // namespace math
} // namespace nla3d

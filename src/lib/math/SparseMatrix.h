// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once

#include "sys.h"
#include "math/Vec.h"

namespace nla3d {
namespace math {

// NOTE: Sparse Matrices entries have indexes started from 1 !

class SparseMatrix;
class SparseSymMatrix;

// TODO: Sparse Matrices in nla3d directly used in MKLs' PARDISO equation solver. That's why we use
// uint32 for value indexes. More general way is to use templates with arbitrary type for entities
// and for indexes.

// Class to hold columns and iofeir (Index Of First Element In the Row) of 3-arrays CSR storage
// format. SparsityInfo can be shared between sparse matrices.
class SparsityInfo {
  public:
    SparsityInfo();
    SparsityInfo(uint32 _nrows, uint32 _ncols, uint32 _max_in_row);
    ~SparsityInfo();

    void reinit(uint32 _nrows, uint32 _ncols, uint32 _max_in_row);
    // add information than (_row, _column) entries has non-zero value
    // _i, _j indexes are started from 1
    void addEntry(uint32 _i, uint32 _j);
    // perform compression procedure. After that positions of non-zero entries can't be changed
    void compress();
    bool isCompressed();

    // return number of entries in the _row (_row > 0)
    uint32 nElementsInRow(uint32 _row);
    // get index in values array (see SparseMatrix implementation) for entry position _i, _j
    uint32 getIndex(uint32 _i, uint32 _j);

    friend class BaseSparseMatrix;
    friend class SparseMatrix;
    friend class SparseSymMatrix;
    friend void matBVprod(SparseSymMatrix &B, const dVec &V, const double coef, dVec &R);
    friend void matBVprod(SparseMatrix &B, const dVec &V, const double coef, dVec &R);
    friend void matBTVprod(SparseMatrix &B, const dVec &V, const double coef, dVec &R);

  private:
    void clear();

    // Data arrays to implement compressed sparse row format with 3 arrays (3-array CSR).
    // Format was implemented by using MKL manual.
    // Here we don't need to know exactly values. We need to store only info about non-zeros
    // entries in matrix.
    // Columns array stores column number of entry
    uint32* columns = nullptr;
    // index of first entry in row: iofeir[_row-1] first entry in _row,
    // iofeir[_row]-1 - last entry in _row
    uint32* iofeir = nullptr;
    // number of non-zero entries in the matrix. Also, the size of values array.
    uint32 numberOfValues = 0;
    uint32 maxInRow = 0;
    bool compressed = false;

    // size of the matrix
    uint32 nRows = 0;
    uint32 nColumns = 0;

    static const uint32 invalid;
};


// Base class for Sparse Matrix in a row-oriented data format. This class doesn't have particular
// meaning. See SparseMatrix and SparseSymMatrix for practical usage.
// The implementation of compressed sparse row format with 3 arrays (3-array CSR).
class BaseSparseMatrix {
  public:
    BaseSparseMatrix();
    BaseSparseMatrix(uint32 _nrows, uint32 _ncolumns, uint32 _max_in_row = 100);
    BaseSparseMatrix(std::shared_ptr<SparsityInfo> spar_info);
    ~BaseSparseMatrix();
    
    // perform compression procedure, after this new entries can't be added to matrix
    void compress();

    // for debug purpose
    void printInternalData (std::ostream& out);

    // zero all entries
    void zero();

    // getters
    double* getValuesArray();
    uint32* getColumnsArray();
    uint32* getIofeirArray();
    uint32 nValues();
    uint32 nRows();
    uint32 nColumns();
    std::shared_ptr<SparsityInfo> getSparsityInfo();
    void setSparsityInfo(std::shared_ptr<SparsityInfo> spar_info);

    bool isCompressed();

  protected:
    // particular values array. Part of compressed sparse format with 3 arrays (3-array CSR).
    double *values = nullptr;

    std::shared_ptr<SparsityInfo> si;
};



class SparseMatrix : public BaseSparseMatrix {
  public:
    SparseMatrix();
    SparseMatrix(uint32 nrows, uint32 ncolumns, uint32 max_in_row = 100);
    SparseMatrix(std::shared_ptr<SparsityInfo> spar_info);
    
    void reinit(uint32 _nrows, uint32 _ncols, uint32 _max_in_row = 100);

    // add non-zero entry to sparse matrix. This should be called before compress(). 
    void addEntry(uint32 _i, uint32 _j);

    // add value to the _i, _j entry. This should be called after compress().
    void addValue(uint32 _i, uint32 _j, double value);

    // debug output
    void print (std::ostream& out);

    // return reference to the entry, if the entry doesn't exist in the matrix - fatal error
    double& operator()(uint32 _i, uint32 _j);

    // return the value of the entry, if the entry doesn't exists - return 0.0
    double value(uint32 _i, uint32 _j) const;

    friend void matBVprod(SparseMatrix &B, const dVec &V, const double coef, dVec &R);
    friend void matBTVprod(SparseMatrix &B, const dVec &V, const double coef, dVec &R);
};


class SparseSymMatrix : public BaseSparseMatrix {
  public:
    SparseSymMatrix();
    SparseSymMatrix(uint32 nrows, uint32 max_in_row = 100);
    SparseSymMatrix(std::shared_ptr<SparsityInfo> spar_info);

    void reinit(uint32 _nrows, uint32 _max_in_row = 100);
    
    // add non-zero entry to sparse matrix. This should be called before compress(). 
    void addEntry(uint32 _i, uint32 _j);

    // add value to the _i, _j entry. This should be called after compress().
    void addValue(uint32 _i, uint32 _j, double value);

    // debug methods
    void print (std::ostream& out);

    // return reference to the entry, if the entry doesn't exist in the matrix - fatal error
    double& operator()(uint32 _i, uint32 _j);

    // return the value of the entry, if the entry doesn't exists - return 0.0
    double value(uint32 _i, uint32 _j) const;

    friend void matBVprod(SparseSymMatrix &B, const dVec &V, const double coef, dVec &R);
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

inline uint32 BaseSparseMatrix::nValues() {
  assert(si);
  return si->numberOfValues;
}

inline uint32 BaseSparseMatrix::nRows() {
  assert(si);
  return si->nRows;
}

inline uint32 BaseSparseMatrix::nColumns() {
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


inline void SparseMatrix::addEntry(uint32 _i, uint32 _j) {
  assert(si);
  si->addEntry(_i, _j);
}

inline void SparseMatrix::addValue(uint32 _i, uint32 _j, double value) {
  this->operator()(_i, _j) += value;
}


inline void SparseSymMatrix::addEntry(uint32 _i, uint32 _j) {
  assert(si);

  // ensure that we work in upper triangle
	if (_i > _j) std::swap(_i, _j);
  si->addEntry(_i, _j);
}

inline void SparseSymMatrix::addValue(uint32 _i, uint32 _j, double value) {
  this->operator()(_i, _j) += value;
}

 
} // namespace math
} // namespace nla3d

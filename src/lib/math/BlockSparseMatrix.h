// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once

#include "sys.h"
#include "math/SparseMatrix.h"

namespace nla3d {
namespace math {


// BlockSparseSymMatrix it is useful wrap over several SparseMatrix (and SparseSymMatrix) when your
// global sparse matrix is represented by different blocks (`nb` is number of blocks). Here is a
// illustration:
// Here is global symmetric matrix A `n` by `n`, it consists of 3 distinct sparse matrix blocks:
//     | A1  | A12 |    
// A = |-----------|,
//     |A12^T| A2  |
//
// where A1 is `m` by `m` symmetric sparse matrix, A2 is `l` by `l` symmetric sparse matrix and A12
// `m` by `l` sparse matrix. `n` = `m` + `l`.
//
// Such matrix consisted of distinct blocks can be created with:
// BlockSparseSymMatrix<2> A(m, n);
//
// Different blocks can be accessed by block() function:
// SparseSymMatrix* A1 = A.block(1);
// SparseSymMatrix* A2 = A.block(2);
// SparseMatrix* A12 = A.block(1, 2);
//
// All others methods like zero(), compress() and others directly trigger the same methods of blocks
// matrices.
//
// Methods addEntry(), addValue() works with global indices (in our example it's from 1 to `n`).
// Under the hood these methods decide in which block translate your request, then modify indices to
// local and call block's method.
template <uint16 nb>
class BlockSparseSymMatrix {
  public:
    BlockSparseSymMatrix(std::initializer_list<uint32> _rows_in_block, uint32 max_in_row = 100);
    BlockSparseSymMatrix(BlockSparseSymMatrix* ex);
    
    // add non-zero entry to sparse matrix. This should be called before compress(). 
    void addEntry(uint32 _i, uint32 _j);

    // add value to the _i, _j entry. This should be called after compress().
    void addValue(uint32 _i, uint32 _j, double value);

    SparseSymMatrix* block(uint16 _i);
    SparseMatrix* block(uint16 _i, uint16 _j);

    void zero();
    void compress();
    bool isCompressed();

    uint32 nRows();

  private:
    void getBlockAndPosition(uint32 _i, uint16* block, uint32* pos);

    SparseSymMatrix diag[nb];
    SparseMatrix upper[nb * (nb + 1) / 2 - nb];

    std::vector<uint32> rows_in_block;
    uint32 total_rows = 0;
    bool compressed = false;
};


template<uint16 nb>
BlockSparseSymMatrix<nb>::BlockSparseSymMatrix(std::initializer_list<uint32> _rows_in_block, uint32 _mat_in_row) {
  rows_in_block = std::vector<uint32>(_rows_in_block);

  assert(rows_in_block.size() == nb);

  for (uint16 i = 0; i < nb; i++) {
    block(i+1)->reinit(rows_in_block[i], _mat_in_row);
    for (uint16 j = i+1; j < nb; j++) {
      block(i+1, j+1)->reinit(rows_in_block[i], rows_in_block[j], _mat_in_row);
    }
    total_rows += rows_in_block[i];
  }
}


template<uint16 nb>
BlockSparseSymMatrix<nb>::BlockSparseSymMatrix(BlockSparseSymMatrix<nb>* ex) {
  // assign sparsity to all underlying matrices
  rows_in_block = ex->rows_in_block; // TODO: need copy here..
  assert(rows_in_block.size() == nb);

  total_rows = ex->total_rows;
  compressed = ex->compressed;

  for (uint16 i = 0; i < nb; i++) {
    block(i+1)->setSparsityInfo(ex->block(i+1)->getSparsityInfo());
    for (uint16 j = i+1; j < nb; j++) {
      block(i+1, j+1)->setSparsityInfo(ex->block(i+1, j+1)->getSparsityInfo());
    }
    total_rows += rows_in_block[i];
  }
}


template<uint16 nb>
inline SparseSymMatrix* BlockSparseSymMatrix<nb>::block(uint16 _i) {
  assert(_i > 0 && _i <= nb);
  return &diag[_i-1];
}


template<uint16 nb>
inline SparseMatrix* BlockSparseSymMatrix<nb>::block(uint16 _i, uint16 _j) {
  CHECK(_i != _j);
  if (_i > _j) std::swap(_i, _j);

  uint16 ind = nb*(_i-1) - (_i-2)*(_i-1)/2 + (_j-_i) - _i;
  assert(ind >= 0 && ind < nb * (nb + 1) / 2 - nb);

  return &upper[ind];
}


template<uint16 nb>
inline void BlockSparseSymMatrix<nb>::getBlockAndPosition(uint32 _i, uint16* block, uint32* pos) {
  assert(_i > 0 && _i <= total_rows);
  *block = 1;
  *pos = _i;
  for (auto& blocknum : rows_in_block) {
    if (*pos <= blocknum) break;
    *pos -= blocknum;
    *block += 1;
  }
}


template<uint16 nb>
void BlockSparseSymMatrix<nb>::addEntry(uint32 _i, uint32 _j) {
  uint16 block_i, block_j;
  uint32 pos_i, pos_j;
  if (_i > _j) std::swap(_i, _j);
  getBlockAndPosition(_i, &block_i, &pos_i);
  getBlockAndPosition(_j, &block_j, &pos_j);

  if (block_i == block_j) {
    block(block_i)->addEntry(pos_i, pos_j);
  } else {
    block(block_i, block_j)->addEntry(pos_i, pos_j);
  }
}


template<uint16 nb>
void BlockSparseSymMatrix<nb>::addValue(uint32 _i, uint32 _j, double value) {
  uint16 block_i, block_j;
  uint32 pos_i, pos_j;
  if (_i > _j) std::swap(_i, _j);
  getBlockAndPosition(_i, &block_i, &pos_i);
  getBlockAndPosition(_j, &block_j, &pos_j);

  if (block_i == block_j) {
    block(block_i)->addValue(pos_i, pos_j, value);
  } else {
    block(block_i, block_j)->addValue(pos_i, pos_j, value);
  }
}

template<uint16 nb>
void BlockSparseSymMatrix<nb>::zero() {
  for (uint16 i = 0; i < nb; i++) {
    diag[i].zero();
  }

  for (uint16 i = 0; i < nb * (nb + 1) / 2 - nb; i++) {
    upper[i].zero();
  }
}


template<uint16 nb>
void BlockSparseSymMatrix<nb>::compress() {
  // NOTE: now this function can be called multiply times. We rely on underlying
  // SparseMatrix::compress()
  for (uint16 i = 0; i < nb; i++) {
    diag[i].compress();
  }

  for (uint16 i = 0; i < nb * (nb + 1) / 2 - nb; i++) {
    upper[i].compress();
  }

  compressed = true;
}

template<uint16 nb>
inline bool BlockSparseSymMatrix<nb>::isCompressed(){
  return compressed;
}

template<uint16 nb>
uint32 BlockSparseSymMatrix<nb>::nRows() {
  return total_rows;
}


} // math
} // nla3d

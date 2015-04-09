// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once

#include "sys.h"
#include <map>

namespace nla3d {
namespace math {

// TODO:
// SparseMatrix in nla3d is primarly used to provide sparse matrix data
// to MKL's PARDISO solver. That's the reason to use uint32 type for indexing.
// It's more general to use size_t for indexes. To manage with it there are several ways:
// 1. Use templated calsses
// 2. Use typedef int SparseIndex_t
// 3. Use #define SparseIndex_t int

struct ValueInfoForSparseMatrix {
	ValueInfoForSparseMatrix() : value(0.0), row(0), column(0) { }
	ValueInfoForSparseMatrix(double v, uint32 r, uint32 c) : value(v), row(r), column(c) { }
	double value;
	uint32 row;
	uint32 column;
};


// general Sparse Matrix in a row-oriented data format
// The implementation of compressed sparse row format with 3 arrays (3-array CSR).
class SparseMatrix {
public:
	SparseMatrix (uint32 nrows, uint32 ncolumns);
	~SparseMatrix ();
  
	void startTraining ();
	void stopTraining (bool copyValuesToMatrix = true);
  bool isInTrainingMode();

	void addValue (uint32 row, uint32 column, double value);

  // TODO: It would be better to use a BLAS routine for sparse matrices for speedup
	double mult_vec_i (double *vec, uint32 n);
	double transpose_mult_vec_i (double *vec, uint32 n);
	void transpose_mult_vec (double *vec, double *res);

	void clear();

  // debug methods
  void printInternalData (std::ostream& out);
  void print (std::ostream& out);

	void zero();

  // getters
	double& operator() (uint32 row, uint32 column);
	double* getValuesArray();
	uint32* getColumnsArray();
	uint32* getIofeirArray();
	uint32 getNumberOfValues();
	uint32 getNumberOfRows();
	uint32 getNumberOfColumns();

private:
  // stuff for training of a sparse matrix
	uint64 getGeneralIndex(uint32 row, uint32 column);
	bool is_training;
  // the training procedure works by mean of std::map container
	std::map<uint64, ValueInfoForSparseMatrix> training_data;
  // counter for number of non-zero element in the every row
	uint32 *numberOfElementsInRow;

  // get index in the values array for [row, column] element
	bool searchIndex (uint32 row, uint32 column, uint32* ind);

  // just a variable to return a pointer to a zero value
	static double dummy;

  // Data arrays to implement compressed sparse row format with 3 arrays (3-array CSR).
  // Format was implemented by using MKL manual.
	double *values;
	uint32 *columns;
  // index of first element in row
	uint32 *iofeir;
  // number of non-zero elements in the matrix. Also, size of values array.
  // If numberOfValues > 0 that means that the matrix is ready for work (a training has been completed already)
	uint32 numberOfValues;

  // size of the matrix
	uint32 numberOfRows;
	uint32 numberOfColumns;
};



// SparseSymmetricMatrix (UpperTriangle)
// The implementation of compressed sparse row format with 3 arrays (3-array CSR).
// TODO: make one template class for SparseMatrix and SparseSymmetricMatrix
class SparseSymmetricMatrix {
public:
	SparseSymmetricMatrix (uint32 nrows);
	~SparseSymmetricMatrix ();

	void startTraining ();
	void stopTraining (bool copyValuesToMatrix = true);
  bool isInTrainingMode();

	void addValue (uint32 row, uint32 column, double value);

  // TODO: It would be better to use a BLAS routine for sparse matrices for speedup
	double mult_vec_i(double *vec, uint32 n);
	void clear();

	void zero();
	void zeroBlock(uint32 n);

  // getters
  double& operator() (uint32 row, uint32 column);
	double* getValuesArray();
	uint32* getColumnsArray();
	uint32* getIofeirArray();
	uint32 getNumberOfValues();
	uint32 getNumberOfRows();

  // debug methods
  void printInternalData (std::ostream& out);
  void print (std::ostream& out);

private:
  // stuff for training of a sparse matrix
	uint64 getGeneralIndex(uint32 row, uint32 column);
	bool is_training;
  // the training procedure works by mean of std::map container
  std::map<uint64, ValueInfoForSparseMatrix> training_data;
  // counter for number of non-zero element in the every row
	uint32 *numberOfElementsInRow;

  // just a variable to return a pointer to a zero value
	static double dummy;

  // get index in the values array for [row, column] element
	bool searchIndex(uint32 row, uint32 column, uint32* ind);

  // Data arrays to implement compressed sparse row format with 3 arrays (3-array CSR).
  // Format was implemented by using MKL manual.
	double* values;
	uint32* columns;
  // index of first element in row
	uint32* iofeir;
  // number of non-zero elements in the matrix. Also, size of values array.
  // If numberOfValues > 0 that means that the matrix is ready for work (a training has been completed already)
	uint32 numberOfValues;
  // size of the matrix
	uint32 numberOfRows;
};

} // namespace math
} // namespace nla3d

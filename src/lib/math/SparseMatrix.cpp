// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "SparseMatrix.h"

namespace nla3d {
namespace math {


const uint32 SparsityInfo::invalid = 0xFFFFFFFF;


SparsityInfo::SparsityInfo(uint32 _nrows, uint32 _ncols, uint32 _max_in_row) {
  assert (_max_in_row > 0);
  nRows = _nrows;
  nColumns = _ncols;


  maxInRow = _max_in_row < nColumns  ? _max_in_row : nColumns;
	columns = new uint32[nRows * maxInRow];
  std::fill_n(columns, nRows * maxInRow, invalid);

	iofeir = new uint32[nRows+1];
  iofeir[0] = 1;
  for (uint32 i = 1; i <= nRows; i++) {
    iofeir[i] = iofeir[i-1] + maxInRow;
  }
}

// row and column positions are started from 1
void SparsityInfo::add(uint32 _i, uint32 _j) {
  assert(iofeir != nullptr);
  assert(columns != nullptr);
  assert(_i > 0 && _i <= nRows);
  assert(_j > 0 && _j <= nColumns);
  assert(compressed == false);
  // NOTE: before compressions columns are not sorted
  uint32 st = iofeir[_i-1] - 1;
  uint32 en = iofeir[_i] - 1;
  for (uint32 i = st; i < en; i++) {
    // if element already in matrix
    if (columns[i] == _j) return;
    if (columns[i] == invalid) {
      columns[i] = _j;
      numberOfValues++;
      return;
    }
  }
  LOG(FATAL) << "maxInRow overflow!";
}


void SparsityInfo::compress() {
  if (compressed) return;
  // count number of added elements
  // uint32 nNonZero = nRows * maxInRow - 
  //                  std::count(&columns[0], &columns[iofeir[nRows]], invalid);
  uint32 next = 0;
  uint32 nextRow = 0;
  uint32* old_columns = columns;
  columns = new uint32[numberOfValues];

  for (uint32 i = 0; i < nRows; i++) {
    for (uint32 j = iofeir[i] - 1; j < iofeir[i+1] - 1; j++) {
      if (old_columns[j] == invalid) break;
      columns[next++] = old_columns[j];
    }
    if ((next - nextRow) > 1) {
      std::sort(&columns[nextRow], &columns[next]);
    }

    iofeir[i] = nextRow + 1;
    nextRow = next;
  }
  // check that we used all space in columns
  assert(next == numberOfValues);
  iofeir[nRows] = next + 1;

  delete[] old_columns;
  compressed = true;
}


uint32 SparsityInfo::getIndex(uint32 _i, uint32 _j) {
  assert(columns);
  assert(iofeir);
	assert(_i > 0 && _i <= nRows);
	assert(_j > 0 && _j <= nColumns);

  uint32 ind;
	uint32 st = iofeir[_i-1] - 1;
	uint32 en = iofeir[_i] - 1;

  if (st == en) return invalid;

  en--;

	while(1) {
		if (en - st == 1) {
			if (columns[st] == _j)
        return st;
			if (columns[en] == _j)
        return en;
			return invalid;
		}
		ind = (uint32) ((en+st)*0.5);
		
		if (columns[ind] == _j)
      return ind;

		if (en == st)
      return invalid;

		if (columns[ind] > _j) {
			en = ind;
    } else {
			st = ind;
    }
	}
  LOG(FATAL) << "What i'm doing here..?";
}



void BaseSparseMatrix::printInternalData(std::ostream& out) {
  assert(si);
  assert(values);
  assert(si->compressed);

	out << "values = {";
	for (uint32 i = 0; i < si->numberOfValues; i++) {
		out << values[i] << "\t";
  }
	out << "}" << std::endl;

	out << "columns = {";
	for (uint32 i = 0; i < si->nRows; i++) {
		for (uint32 j = si->iofeir[i] - 1; j < si->iofeir[i + 1] - 1; j++) {
			out << si->columns[j] << "\t";
    }
	}
	out << "}" << std::endl;

	out << "iofeir = {";
	for (uint32 i = 0; i <= si->nRows; i++) {
		out << si->iofeir[i] << "\t";
  }
	out << "}" << std::endl;
}

BaseSparseMatrix::BaseSparseMatrix(uint32 nrows, uint32 ncolumns, uint32 max_in_row) {
  // create it's own SparsityInfo instance
  si = std::shared_ptr<SparsityInfo>(new SparsityInfo(nrows, ncolumns, max_in_row));
}

BaseSparseMatrix::BaseSparseMatrix(std::shared_ptr<SparsityInfo> spar_info) {
  assert(spar_info);
  // use already exist SparsityInfo instance. It can in compressed = false state or already compressed = true state
  si = spar_info;
}

BaseSparseMatrix::~BaseSparseMatrix() {
	clear();
}

void BaseSparseMatrix::compress() {
  assert(si);
  if (si->compressed == true) {
    if (values)
      // compressed has been finished and values allocated - nothing to do..
      return;
  } else {
    si->compress();
  }
  
  if (values) delete[] values;
  values = new double[si->numberOfValues];
  zero();
}

void BaseSparseMatrix::clear() {
	if (values) {
    delete[] values;
    values = nullptr;
  }
}


SparseMatrix::SparseMatrix(uint32 nrows, uint32 ncolumns, uint32 max_in_row) :
    BaseSparseMatrix(nrows, ncolumns, max_in_row)
{

}

SparseMatrix::SparseMatrix(std::shared_ptr<SparsityInfo> spar_info) :
    BaseSparseMatrix(spar_info)
{

}

double SparseMatrix::mult_vec_i(double *vec, uint32 i) {
  assert(si);
  assert(si->compressed);
  assert(values);
  assert(vec);
  assert(i > 0 && i <= si->nRows);

	double res = 0.0;

  //if no elements in the current row return zero res
	if (si->iofeir[i]-si->iofeir[i-1] == 0)
		return 0.0;

	uint32 st = si->iofeir[i-1] - 1;
	uint32 en = si->iofeir[i] - 2;

	for (uint32 j = st; j <= en; j++)
    res += values[j]*vec[si->columns[j]-1];

	return res;
}

double SparseMatrix::transpose_mult_vec_i(double *vec, uint32 i) {
  assert(si);
  assert(si->compressed);
  assert(values);
  assert(vec);
  assert(i > 0 && i <= si->nColumns);

	double res = 0.0;
	const double eps = 1e-20;

	uint32 index;
	for (uint32 j = 0; j < si->nRows; j++) {
		if (fabs(vec[j]) > eps) {
      index = si->getIndex(j + 1, i);
			if (index != SparsityInfo::invalid)
				res += values[index] * vec[j];
    }
  }
	return res;
}

void SparseMatrix::transpose_mult_vec(double *vec, double *res) {
  assert(si);
  assert(si->compressed);
  assert(values);
	assert(vec && res);

  std::fill(&res[0], &res[si->nColumns], 0.0);

	for (uint32 i = 0; i < si->nRows; i++) {
		if (si->iofeir[i+1]-si->iofeir[i] == 0) {
			continue;
    }
		uint32 st = si->iofeir[i] - 1;
		uint32 en = si->iofeir[i + 1] - 2;
		for (uint32 j = st; j <= en; j++) {
      res[si->columns[j] - 1] += values[j] * vec[i];
    }
	}
}

void SparseMatrix::print(std::ostream& out) {
  assert(si);
  assert(si->compressed);
  assert(values);

	uint32 ind;
	for (uint32 i = 1; i <= si->nRows; i++) {
		out << "[";
		for (uint32 j = 1; j <= si->nColumns; j++) {
      out << value(i, j) << '\t';
    }
		out << "]" << std::endl;
	}
}

double& SparseMatrix::operator() (uint32 _i, uint32 _j) {
  assert(values);
  assert(si);

	uint32 index = si->getIndex(_i, _j);
  if (index == SparsityInfo::invalid) {
    LOG(FATAL) << "The position(" << _i << ", " << _j << ") is absent in the matrix";
  }
  return values[index];
}

double SparseMatrix::value(uint32 _i, uint32 _j) const {
  assert(values);
  assert(si);

	uint32 index = si->getIndex(_i, _j);
  if (index == SparsityInfo::invalid) {
    return 0.0;
  }
  return values[index];
}


SparseSymMatrix::SparseSymMatrix(uint32 nrows, uint32 max_in_row) :
    BaseSparseMatrix(nrows, nrows, max_in_row)
{
  assert(si);
  // NOTE: need to add diagonal elements because MKL PARDISO need it anyway
  for (uint32 i = 1; i <= si->nRows; i++) {
    si->add(i, i);
  }
}

SparseSymMatrix::SparseSymMatrix(std::shared_ptr<SparsityInfo> spar_info) :
    BaseSparseMatrix(spar_info)
{
  // TODO: we need to check that spar_info meets symmetry requirements
  // Let's at leas check that is is rectangular

  assert(spar_info->nRows == spar_info->nColumns);
}

double SparseSymMatrix::mult_vec_i(double *vec, uint32 i) {
  assert(si);
  assert(si->compressed);
  assert(values);
  assert(vec);
  assert(i > 0 && i <= si->nRows);

	double res = 0.0;
	double eps = 1e-20;
	uint32 index;

  // walk on lower triangle
	for (uint32 j = 1; j <= i; j++)
		if (fabs(vec[j-1]) > eps)
      res += value(j, i) * vec[j - 1];

  // walk on upper triangle
	for (uint32 j = i + 1; j <= si->nColumns; j++)
		if (fabs(vec[j-1]) > eps)
      res += value(i, j) * vec[j - 1];

	return res;
}


void SparseSymMatrix::print(std::ostream& out) {
  assert(si);
  assert(si->compressed);
  assert(values);

	uint32 ind;
	for (uint32 i = 1; i <= si->nRows; i++) {
		// out << "[";
		for (uint32 j = 1; j <= si->nColumns; j++) {
      out << value(i, j) << '\t';
    }
		// out << "]" << std::endl;
    out << std::endl;
	}
}

double& SparseSymMatrix::operator() (uint32 _i, uint32 _j) {
  assert(values);
  assert(si);

  // ensure that we work in upper triangle
	if (_i > _j) std::swap(_i, _j);

	uint32 index = si->getIndex(_i, _j);
  if (index == SparsityInfo::invalid) {
    LOG(FATAL) << "The position(" << _i << ", " << _j << ") is absent in the matrix";
  }
  return values[index];
}


double SparseSymMatrix::value(uint32 _i, uint32 _j) const {
  assert(values);
  assert(si);

  // ensure that we work in upper triangle
	if (_i > _j) std::swap(_i, _j);

	uint32 index = si->getIndex(_i, _j);
  if (index == SparsityInfo::invalid) {
    return 0.0;
  }
  return values[index];
}

  
} // namespace math
} // namespace nla3d

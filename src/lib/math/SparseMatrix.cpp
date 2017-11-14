// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "SparseMatrix.h"

namespace nla3d {
namespace math {


const uint32 SparsityInfo::invalid = 0xFFFFFFFF;

SparsityInfo::SparsityInfo() {

}


SparsityInfo::SparsityInfo(uint32 _nrows, uint32 _ncols, uint32 _max_in_row) {
  reinit(_nrows, _ncols, _max_in_row);
}


SparsityInfo::~SparsityInfo() {
  clear();
}


void SparsityInfo::clear() {
  if (iofeir) {
    delete[] iofeir;
  }
  iofeir = nullptr;

  if (columns) {
    delete[] columns;
  }
  columns = nullptr;

  compressed = false;
  numberOfValues = 0;
}


void SparsityInfo::reinit(uint32 _nrows, uint32 _ncols, uint32 _max_in_row) {
  assert (_max_in_row > 0);
  clear();

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
void SparsityInfo::addEntry(uint32 _i, uint32 _j) {
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
    // have a space to store new entry in the row
    if (columns[i] == invalid) {
      columns[i] = _j;
      numberOfValues++;
      return;
    }
  }
  // no more room for new entry
  LOG(FATAL) << "maxInRow overflow!";
}


void SparsityInfo::compress() {
  if (compressed) return;

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


void BaseSparseMatrix::writeCoordinateTextFormat(std::ostream& out) {
  assert(si);
  assert(values);
  assert(si->compressed);
  // we need to return back old preferences after usage
  auto old_precision = out.precision(15);
  auto old_flags = out.setf(std::ios_base::scientific, std::ios_base::floatfield);
  // write header: nRows, nColumns, numberOfValues
  out << si->nRows << ' ' << si->nColumns << ' ' << si->numberOfValues << std::endl;
  uint32 total = 0;
  for (uint32 i = 1; i <= si->nRows; i++) {
    for (uint32 j = si->iofeir[i - 1] - 1; j < si->iofeir[i] - 1; j++) {
        out << i << ' ' << si->columns[j] << ' ' << values[j] << std::endl;
        total++;
    }
  }
  assert(total == si->numberOfValues);
  // return back old setup
  out.precision(old_precision);
  out.setf(old_flags);
}


void BaseSparseMatrix::readCoordinateTextFormat(std::istream& in) {
  std::vector<SparseEntry> entries;
  uint32 _nrows, _ncols, _nvalues;
  // read header: nRows, nColumns, numberOfValues
  in >> _nrows >> _ncols >> _nvalues;
  entries.reserve(_nvalues);
  for (uint32 i = 1; i <= _nvalues; i++) {
    SparseEntry entry;
    in >> entry.i >> entry.j >> entry.v;
    entries.push_back(entry);
  }
  reinit(_nrows, _ncols, entries);
}


bool BaseSparseMatrix::compare(const BaseSparseMatrix& op2, double th) {
  // this function conduct strict comparison: matrices should have the same sparsity and the same
  // non-zero values
  BaseSparseMatrix& op1 = *this;
  assert(op1.si && op2.si);
  assert(op1.values && op2.values);
  assert(op1.si->compressed && op2.si->compressed);

  // first round. compare shapes
  if (op1.nRows() != op2.nRows() || op1.nColumns() != op2.nColumns()) {
    return false;
  }

  // second round. compare number of values
  if (op1.nValues() != op2.nValues()) {
    return false;
  }

  // third round. compare sparsity patterns
  for (uint32 i = 0; i <= op1.nRows(); i++) {
    if (op1.si->iofeir[i] != op2.si->iofeir[i]) {
      return false;
    }
  }
  for (uint32 i = 0; i < op1.nValues(); i++) {
    if (op1.si->columns[i] != op2.si->columns[i]) {
      return false;
    }
  }

  // fourth round. compare values with thresholds
  for (uint32 i = 0; i < op1.nValues(); i++) {
    if (fabs(op1.values[i] - op2.values[i]) > th) {
      return false;
    }
  }

  // all rounds are passed. The matrices are the same
  return true;
}


BaseSparseMatrix::BaseSparseMatrix() {

}

BaseSparseMatrix::BaseSparseMatrix(uint32 _nrows, uint32 _ncolumns, uint32 _max_in_row) {
  // create it's own SparsityInfo instance
  si = std::shared_ptr<SparsityInfo>(new SparsityInfo(_nrows, _ncolumns, _max_in_row));
}

BaseSparseMatrix::BaseSparseMatrix(std::shared_ptr<SparsityInfo> spar_info) {
  setSparsityInfo(spar_info);
}

BaseSparseMatrix::~BaseSparseMatrix() {
  if (values) {
      delete[] values;
      values = nullptr;
    }
}


void BaseSparseMatrix::reinit(uint32 _nrows, uint32 _ncols, const std::vector<SparseEntry>& entries) {
  // after this procedure we will obtain matrix in already compressed state
  // TODO: this is not optimal way to fill new Sparse Matrix, actually we can allocate exactly
  // number of non-zeros at once, as far as we know all entries before initialization
  // NOTE: For SparseSymMatrix the caller should be sure that all entries are in upper triangle.
  //
  // in case of all zeros matrix
  if (entries.size() == 0) {
    si = std::shared_ptr<SparsityInfo>(new SparsityInfo(_nrows, _ncols, 1));
    compress();
    return;
  }
  // calculate number of entries in every row
  std::vector<uint32> entriesInRow(_nrows, 0);
  for (auto& v : entries) {
    entriesInRow[v.i - 1] += 1;
  }
  uint32 _max_in_row = *std::max_element(entriesInRow.begin(), entriesInRow.end());
  // init new SparsityInfo
  si = std::shared_ptr<SparsityInfo>(new SparsityInfo(_nrows, _ncols, _max_in_row));

  // add non-zero entries into SparsityInfo
  for (auto& v : entries) {
    si->addEntry(v.i, v.j);
  }

  // compress matrix
  compress();

  // write non-zero values into matrix
  for (auto& v : entries) {
    //addValue(v.i, v.j, v.v);
    uint32 index = si->getIndex(v.i, v.j);
    if (index == SparsityInfo::invalid) {
      LOG(FATAL) << "The position(" << v.i << ", " << v.j << ") is absent in the matrix";
    }
    values[index] = v.v;
  }
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


void BaseSparseMatrix::setSparsityInfo(std::shared_ptr<SparsityInfo> spar_info) {
  assert(spar_info);
  // use already exist SparsityInfo instance. It can in compressed = false state or already
  // compressed = true state
  si = spar_info;
  // if we add SparsityInfo with alreade fixed number of entries we are ready to allocate values
  if (si->compressed == true) 
    compress();
}


SparseMatrix::SparseMatrix() {

}


SparseMatrix::SparseMatrix(uint32 _nrows, uint32 _ncolumns, uint32 _max_in_row) :
    BaseSparseMatrix(_nrows, _ncolumns, _max_in_row)
{

}

SparseMatrix::SparseMatrix(std::shared_ptr<SparsityInfo> spar_info) :
    BaseSparseMatrix(spar_info)
{

}


void SparseMatrix::reinit(uint32 _nrows, uint32 _ncols, uint32 _max_in_row) {
  si = std::shared_ptr<SparsityInfo>(new SparsityInfo(_nrows, _ncols, _max_in_row));
  if (values) {
    delete[] values;
    values = nullptr;
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


double& SparseMatrix::operator()(uint32 _i, uint32 _j) {
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


SparseSymMatrix::SparseSymMatrix() {

}


SparseSymMatrix::SparseSymMatrix(uint32 _nrows, uint32 _max_in_row) :
    BaseSparseMatrix(_nrows, _nrows, _max_in_row)
{
  assert(si);
  // NOTE: need to add diagonal elements because MKL PARDISO need it anyway
  for (uint32 i = 1; i <= si->nRows; i++) {
    si->addEntry(i, i);
  }
}

SparseSymMatrix::SparseSymMatrix(std::shared_ptr<SparsityInfo> spar_info) :
    BaseSparseMatrix(spar_info)
{
  // TODO: we need to check that spar_info meets symmetry requirements
  // Let's at least check that is is rectangular

  assert(spar_info->nRows == spar_info->nColumns);
}

void SparseSymMatrix::reinit(uint32 _nrows, uint32 _max_in_row) {
  si = std::shared_ptr<SparsityInfo>(new SparsityInfo(_nrows, _nrows, _max_in_row));
  if (values) {
    delete[] values;
    values = nullptr;
  }
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


void matBVprod(SparseSymMatrix &B, const dVec &V, const double coef, dVec &R) {
  assert(B.si);
  assert(B.si->compressed);
  assert(B.values);
  assert(B.nRows() == V.size());
  assert(R.size() >= B.nRows());
  // TODO: Try to use BLAS routines and measure speedup

	const double eps = 1e-20;

  for (uint32 i = 1; i <= B.nRows(); i++) {
    // walk on lower triangle
    for (uint32 j = 1; j <= i; j++)
      if (fabs(V[j-1]) > eps)
        R[i-1] += B.value(j, i) * V[j-1] * coef;

    // walk on upper triangle
    // TODO: could be done in more efficient way
    for (uint32 j = i + 1; j <= B.nRows(); j++)
      if (fabs(V[j-1]) > eps)
        R[i-1] += B.value(i, j) * V[j-1] * coef;
  }
}


void matBVprod(SparseMatrix &B, const dVec &V, const double coef, dVec &R) {
  assert(B.si);
  assert(B.si->compressed);
  assert(B.values);
  assert(B.nColumns() == V.size());
  assert(R.size() >= B.nRows());
  // TODO: Try to use BLAS routines and measure speedup

  for (uint32 i = 1; i <= B.nRows(); i++) {

    //if no elements in the current row return zero res
    if (B.si->iofeir[i] - B.si->iofeir[i-1] == 0)
      continue;

    uint32 st = B.si->iofeir[i-1] - 1;
    uint32 en = B.si->iofeir[i] - 2;

    for (uint32 j = st; j <= en; j++)
      R[i-1] += B.values[j] * V[B.si->columns[j]-1] * coef;
  }
}

void matBTVprod(SparseMatrix &B, const dVec &V, const double coef, dVec &R) {
  assert(B.si);
  assert(B.si->compressed);
  assert(B.values);
  assert(B.nRows() == V.size());
  assert(R.size() >= B.nColumns());
  // TODO: Try to use BLAS routines and measure speedup

	for (uint32 i = 1; i <= B.si->nRows; i++) {
		if (B.si->iofeir[i] - B.si->iofeir[i-1] == 0) {
			continue;
    }
		uint32 st = B.si->iofeir[i-1] - 1;
		uint32 en = B.si->iofeir[i] - 2;
		for (uint32 j = st; j <= en; j++)
      R[B.si->columns[j] - 1] += B.values[j] * V[i-1] * coef;
	}
}
  
} // namespace math
} // namespace nla3d

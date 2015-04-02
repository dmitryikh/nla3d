// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "SparseMatrix.h"

namespace nla3d {
namespace math {

double SparseMatrix::dummy = 0.0;
double SparseSymmetricMatrix::dummy = 0.0;


SparseMatrix::SparseMatrix(uint32 nrows, uint32 ncolumns) : numberOfRows(nrows), numberOfColumns(ncolumns) {
	is_training = false;
	values = NULL;
	columns = NULL;
	numberOfElementsInRow = NULL;
	numberOfValues = 0;
}

SparseMatrix::~SparseMatrix() {
	clear();
}

void SparseMatrix::startTraining() {
	assert(numberOfRows && numberOfColumns);
	clear();
	if (is_training) {
		warning("SparseMatrix::startTraining: strart trainig before ending the \
        previus one. The previus one has been deleted.");
	}
	is_training = true;
	numberOfElementsInRow = new uint32[numberOfRows];
	memset(numberOfElementsInRow,0,sizeof(uint32)*numberOfRows);
	training_data.clear();
}

void SparseMatrix::addValue(uint32 row, uint32 column, double value) {
	assert(row > 0 && row < numberOfRows+1);
	assert(column > 0 && column < numberOfColumns+1);
	uint64 gi = getGeneralIndex(row, column);
	if (is_training) {
		if (training_data.count(gi) == 0) {
			training_data[gi] = ValueInfoForSparseMatrix(value, row, column);
			numberOfElementsInRow[row-1]++;
		}	else {
			training_data[gi].value += value;
    }
	} else {
		uint32 index;
		if (searchIndex(row, column, &index)) {
			values[index] += value;
    } else {
			warning("SparseMatrix::addValue: the position(%d,%d) is absent on current \
          matrix", row, column);
    }
	}
}

void SparseMatrix::stopTraining(bool copyValuesToMatrix) {
  //TODO: need a check for overflow
	numberOfValues = static_cast<uint32> (training_data.size());
	values = new double[numberOfValues];
	columns = new uint32[numberOfValues];
	iofeir = new uint32[numberOfRows+1];
	iofeir[0] = 1;
	for (uint32 i=0; i < numberOfRows; i++) {
		iofeir[i+1] = iofeir[i]+numberOfElementsInRow[i];
	}

	std::map<uint64, ValueInfoForSparseMatrix>::iterator p = training_data.begin();
	uint32 counter = 0;
	while (p!=training_data.end()) {
		if (copyValuesToMatrix) {
			values[counter] = p->second.value;
    } else {
			values[counter] = 0.0;
    }
		columns[counter] = p->second.column;
		p++;
		counter++;
	}
	assert(counter == numberOfValues);
	// delete all training information (not need any more)
	if (numberOfElementsInRow) {
    delete[] numberOfElementsInRow;
    numberOfElementsInRow = NULL;
  }
	training_data.clear();
	is_training = false;
}

bool SparseMatrix::isInTrainingMode() {
  return is_training;
}

void SparseMatrix::clear() {
	numberOfValues = 0;
	if (values) {
    delete[] values;
    values = NULL;
  }
	if (columns) {
    delete[] columns;
    columns = NULL;
  }
	if (numberOfElementsInRow) {
    delete[] numberOfElementsInRow;
    numberOfElementsInRow = NULL;
  }
}

void SparseMatrix::zero() {
  //TODO: do we need to work with zero sized matrices?
	//assert(numberOfValues);
	memset(values, 0, sizeof(double)*numberOfValues);
}

uint64 SparseMatrix::getGeneralIndex(uint32 row, uint32 column) {
	assert(numberOfColumns);
	return ((uint64)row - 1) * (uint64) numberOfColumns + (uint64)column - (uint64)1;
}

bool SparseMatrix::searchIndex(uint32 row, uint32 column, uint32* ind) {
	assert(numberOfValues);
	assert(row > 0 && row < numberOfRows+1);
	assert(column > 0 && column < numberOfColumns+1);

	if (iofeir[row]-iofeir[row-1] == 0) {
		return false;
  }

	uint32 st = iofeir[row-1]-1;
	uint32 en = iofeir[row]-2;
	while(1) {
		if (en - st == 1) {
			if (columns[st] == column)  {
				*ind = st;
				break;
			}
			if (columns[en] == column) {
				*ind = en;
				break;
			}
			return false;
		}
		*ind = (uint32) ((en+st)*0.5);
		
		if (columns[*ind] == column) {
      break;
    }
		if (en == st) {
			return false;
    }

		if (columns[*ind] > column) {
			en = *ind;
    } else {
			st = *ind;
    }
	}
	return true;
}

double& SparseMatrix::operator() (uint32 row, uint32 column) {
	uint32 index;
	if (searchIndex(row, column, &index)) {
		return values[index];
  } else {
		warning("SparseMatrix::operator(): the position(%d,%d) is absent on current \
        matrix", row, column);
    return dummy;
  }
}


void SparseMatrix::printInternalData(std::ostream& out) {
	out << "values = {";
	for (uint32 i=0; i < numberOfValues; i++) {
		out << values[i] << "\t";
  }
	out << "}" << std::endl;

	out << "columns = {";
	for (uint32 i=0; i < numberOfRows; i++) {
		for (uint32 j=iofeir[i]; j < iofeir[i+1]; j++) {
			out << columns[j-1] << "\t";
    }
	}
	out << "}" << std::endl;

	out << "iofeir = {";
	for (uint32 i=0; i < numberOfRows+1; i++) {
		out << iofeir[i] << "\t";
  }
	out << "}" << std::endl;
}

void SparseMatrix::print(std::ostream& out) {
	uint32 ind;
	for (uint32 i = 0; i < numberOfRows; i++) {
		out << "[";
		for (uint32 j = 0; j < numberOfColumns; j++) {
			if (searchIndex(i+1, j+1, &ind)) {
				out << values[ind] << "\t";
      } else {
				out << "0\t";
      }
    }
		out << "]" << std::endl;
	}
}

double SparseMatrix::mult_vec_i(double *vec, uint32 n) {
	double res = 0.0;
	double eps = 1e-20;

	if (iofeir[n]-iofeir[n-1] == 0) {
    //no elements in the current row
		return 0.0;
  }
	uint32 st = iofeir[n-1]-1;
	uint32 en = iofeir[n]-2;

	for (uint32 j = st; j <= en; j++) {
    res += values[j]*vec[columns[j]-1];
  }
	return res;
}

double SparseMatrix::transpose_mult_vec_i(double *vec, uint32 n) {
	double res = 0.0;
	double eps = 1e-20;
	uint32 index;
	for (uint32 j = 0; j < numberOfRows; j++) {
		if (fabs(vec[j]) > eps) {
			if (searchIndex(j+1, n, &index)) {
				double ve = vec[j];
				double va = values[index];
				res += values[index]*vec[j];
			}
    }
  }
	return res;
}

void SparseMatrix::transpose_mult_vec(double *vec, double *res) {
	assert(vec && res);

	memset(res,0,sizeof(double)*numberOfColumns);
	double eps = 1e-20;

	for (uint32 i = 0; i < numberOfRows; i++) {
		if (iofeir[i+1]-iofeir[i] == 0) {
			continue;
    }
		uint32 st = iofeir[i]-1;
		uint32 en = iofeir[i+1]-2;
		for (uint32 j=st; j <= en; j++) {
      res[columns[j]-1] += values[j]*vec[i];
    }
	}
}

double* SparseMatrix::getValuesArray() {
  assert(numberOfValues);
  return values;
}

uint32* SparseMatrix::getColumnsArray() {
  assert(numberOfValues);
  return columns;
}

uint32* SparseMatrix::getIofeirArray() {
  assert(numberOfValues);
  return iofeir;
}

uint32 SparseMatrix::getNumberOfValues() {
  return numberOfValues;
}

uint32 SparseMatrix::getNumberOfRows() {
  return numberOfRows;
}

uint32 SparseMatrix::getNumberOfColumns() {
  return numberOfColumns;
}




SparseSymmetricMatrix::SparseSymmetricMatrix(uint32 nrows) : numberOfRows(nrows) {
	is_training = false;
	values = NULL;
	columns = NULL;
	numberOfElementsInRow = NULL;
	numberOfValues = 0;
}

SparseSymmetricMatrix::~SparseSymmetricMatrix() {
	clear();
}

void SparseSymmetricMatrix::startTraining() {
	assert(numberOfRows);
	clear();
	if (is_training) {
		warning("SparseSymmetricMatrix::startTraining: strart trainig before ending \
       the previus one. The previus one has been deleted.");
	}
	is_training = true;
	if (numberOfElementsInRow) {
    delete[] numberOfElementsInRow;
  }

	numberOfElementsInRow = new uint32[numberOfRows];
	memset(numberOfElementsInRow,0,sizeof(uint32)*numberOfRows);
	training_data.clear();
  // TODO: not a general thing: we firstly add all diagonal elements
  // It's needed for PARDISO solver
	for (uint32 i=0;i<numberOfRows;i++) {
		addValue(i+1, i+1, 0.0);
	}
}

void SparseSymmetricMatrix::addValue(uint32 row, uint32 column, double value) {
	assert(row > 0 && row < numberOfRows+1);
	assert(column > 0 && column < numberOfRows+1);
	if (row > column) {
    std::swap(column,row);
  }
	uint64 gi = getGeneralIndex(row, column);
	if (is_training) {
		if (training_data.count(gi) == 0) {
			training_data[gi] = ValueInfoForSparseMatrix(value, row, column);
			numberOfElementsInRow[row-1]++;
		} else {
			training_data[gi].value += value;
    }
	} else {
		uint32 index;
		if (searchIndex(row, column, &index)) {
			values[index] += value;
    } else {
			warning("SparseSymmetricMatrix::addValue: the position(%d,%d) is absent on \
          current matrix", row, column);
    }
	}
}

void SparseSymmetricMatrix::stopTraining(bool copyValuesToMatrix) {
  //TODO: need a check for overflow
	numberOfValues = static_cast<uint32> (training_data.size());
	values = new double[numberOfValues];
	columns = new uint32[numberOfValues];
	iofeir = new uint32[numberOfRows+1];
	iofeir[0] = 1;
	for (uint32 i = 0; i < numberOfRows; i++) {
		iofeir[i+1] = iofeir[i]+numberOfElementsInRow[i];
	}

	std::map<uint64, ValueInfoForSparseMatrix>::iterator p = training_data.begin();
	uint32 counter = 0;
	while (p != training_data.end()) {
		if (copyValuesToMatrix) {
			values[counter] = p->second.value;
    } else {
			values[counter] = 0.0;
    }
		columns[counter] = p->second.column;
		p++;
		counter++;
	}
	assert(counter == numberOfValues);
	// delete all training inforamtion (Don't need more)
	if (numberOfElementsInRow) {
    delete[] numberOfElementsInRow;
    numberOfElementsInRow = NULL;
  }
	training_data.clear();
	is_training = false;
}

bool SparseSymmetricMatrix::isInTrainingMode() {
  return is_training;
}

void SparseSymmetricMatrix::clear() {
	numberOfValues = 0;
	if (values) {
    delete[] values;
    values = NULL;
  }
	if (columns) {
    delete[] columns;
    columns = NULL;
  }
	if (numberOfElementsInRow) {
    delete[] numberOfElementsInRow;
    numberOfElementsInRow = NULL;
  }
}

void SparseSymmetricMatrix::zero() {
	assert(numberOfValues);
	memset(values, 0, sizeof(double)*numberOfValues);
}

uint64 SparseSymmetricMatrix::getGeneralIndex(uint32 row, uint32 column) {
	assert(numberOfRows);
  return (uint64)numberOfRows * ((uint64)row-1)-((uint64)row-2)*((uint64)row-1)/2+1+(uint64)column-(uint64)row;
}

bool SparseSymmetricMatrix::searchIndex(uint32 row, uint32 column, uint32* ind) {
	assert(numberOfValues);
	assert(row > 0 && row < numberOfRows+1);
	assert(column > 0 && column < numberOfRows+1);
	if (row > column) {
    std::swap(column,row);
  }

	uint32 st = iofeir[row-1]-1;
	uint32 en = iofeir[row]-2;
	while(1) {
		if (en - st == 1) {
			if (columns[st] == column)  {
				*ind = st;
				break;
			}
			if (columns[en] == column) {
				*ind = en;
				break;
			}
			return false;
		}
		*ind = (uint32) ((en+st)*0.5);
		
		if (columns[*ind] == column) {
      break;
    }
		if (en == st) {
			return false;
    }

		if (columns[*ind] > column) {
			en = *ind;
    } else {
			st = *ind;
    }
	}
	return true;
}

double& SparseSymmetricMatrix::operator() (uint32 row, uint32 column) {
	uint32 index;
	if (searchIndex(row, column, &index)) {
		return values[index];
  } else {
		warning("SparseSymmetricMatrix::operator(): the position(%d,%d) is absent \
        on current matrix", row, column);
    return dummy;
  }
}


void SparseSymmetricMatrix::printInternalData(std::ostream& out) {
	out << "values = {";
	for (uint32 i=0; i < numberOfValues; i++) {
		out << values[i] << "\t";
  }
	out << "}" << std::endl;

	out << "columns = {";
	for (uint32 i=0; i < numberOfRows; i++) {
		for (uint32 j=iofeir[i]; j < iofeir[i+1]; j++) {
			out << columns[j-1] << "\t";
    }
	}
	out << "}" << std::endl;

	out << "iofeir = {";
	for (uint32 i=0; i < numberOfRows+1; i++) {
		out << iofeir[i] << "\t";
  }
	out << "}" << std::endl;
}

void SparseSymmetricMatrix::print (std::ostream& out) {
	uint32 ind;
	for (uint32 i=0; i < numberOfRows; i++) {
		for (uint32 j=0; j < numberOfRows;j++) {
			if (searchIndex(i+1, j+1, &ind)) {
				out << values[ind];
      } else {
				out << "0";
      }
			if (j != numberOfRows-1) {
				out << "\t";
      }
		}
		out << std::endl;
	}
}

double SparseSymmetricMatrix::mult_vec_i(double *vec, uint32 n) {
	double res = 0.0;
	double eps = 1e-20;
	uint32 index;

	for (uint32 j=0; j < numberOfRows; j++) {
		if (fabs(vec[j]) > eps) {
			if (searchIndex(n, j+1, &index)) {
				double ve = vec[j];
				double va = values[index];
				res += values[index]*vec[j];
			}
    }
  }
	return res;
}

void SparseSymmetricMatrix::zeroBlock(uint32 n) {
	assert(n > 0 && n < numberOfRows+1);
	assert(numberOfValues);
	for (uint32 i=0; i < n; i++) {
		uint32 st = iofeir[i]-1;
		uint32 en = iofeir[i+1]-2;
		while (columns[en] > n) {
			en--;
    }
		memset(&values[st],0,sizeof(double)*(en-st+1));
	}
}

double* SparseSymmetricMatrix::getValuesArray() {
  assert(numberOfValues);
  return values;
}

uint32* SparseSymmetricMatrix::getColumnsArray() {
  assert(numberOfValues);
  return columns;
}

uint32* SparseSymmetricMatrix::getIofeirArray() {
  assert(numberOfValues);
  return iofeir;
}

uint32 SparseSymmetricMatrix::getNumberOfValues() {
  return numberOfValues;
}

uint32 SparseSymmetricMatrix::getNumberOfRows()
{
  return numberOfRows;
}

} // namespace math
} // namespace nla3d

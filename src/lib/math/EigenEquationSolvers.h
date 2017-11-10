// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#pragma once

#include "math/SparseMatrix.h"
#include "math/EquationSolver.h"

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

namespace nla3d {
namespace math {


// TODO: LLT-based method currently don't work with RowMajor sparse matrices.
// typedef Eigen::SimplicialLDLT<MatrixType, Eigen::Upper> SolverType;
// typedef Eigen::SimplicialLLT<MatrixType, Eigen::Upper> SolverType;


class EigenCGEquationSolver : public EquationSolver {
public:
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int> MatrixType;
  typedef Eigen::VectorXd VectorType;
  typedef Eigen::VectorXi IndexVectorType;
  typedef Eigen::ConjugateGradient<MatrixType, Eigen::Upper> SolverType;
 
  void solveEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns) override {
    factorizeEquations(matrix);
    substituteEquations(matrix, rhs, unknowns);
  }

  void factorizeEquations(math::SparseSymMatrix* matrix) override {
    // Indexind of elements in math::Sparse*Matrix started from 1 (for MKL's PARDISO input)
    // We store copy of sparse data with indexing for Eigen (started from 0)
    m_iofeir.resize(matrix->nRows() + 1);
    int* const p_iofeir = (int* const) matrix->getIofeirArray();
    for (size_t i = 0; i < m_iofeir.rows(); i++) {
      m_iofeir[i] = p_iofeir[i] - 1;
    }

    m_columns.resize(matrix->nValues());
    int* const p_columns = (int* const) matrix->getColumnsArray();
    for (size_t i = 0; i < m_columns.rows(); i++) {
      m_columns[i] = p_columns[i] - 1;
    }
  };

  void substituteEquations(math::SparseSymMatrix* matrix, double* rhs, double* unknowns) override {
    assert(matrix && rhs && unknowns);
    assert(matrix->nValues() == m_columns.rows());
    assert(matrix->nRows() + 1 == m_iofeir.rows());
    Eigen::Map<const MatrixType> emat( matrix->nRows()
                   , matrix->nColumns()
                   , matrix->nValues()
                   , m_iofeir.data()
                   , m_columns.data()
                   , matrix->getValuesArray()
                   );
    Eigen::Map<VectorType> b(rhs, matrix->nRows());
    Eigen::Map<VectorType> x(unknowns, matrix->nRows());
    SolverType solver;
    x = solver.compute(emat).solve(b);
  };
private:
  IndexVectorType m_iofeir;
  IndexVectorType m_columns;
};

} //math
} //nla3d

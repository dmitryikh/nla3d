// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "math/EquationSolver.h"
#include "math/SparseMatrix.h"
#include <mkl.h>

namespace nla3d {

namespace math {

EquationSolver* defaultEquationSolver = new PARDISO_equationSolver;

void EquationSolver::setSymmetric (bool symmetric) {
  isSymmetric = symmetric;
}

void EquationSolver::setPositive (bool positive) {
  isPositive = positive;
}

PARDISO_equationSolver::~PARDISO_equationSolver () {
  releasePARDISO();
}

void PARDISO_equationSolver::solveEquations (math::SparseSymmetricMatrix* matrix, double* rhs, double* unknowns) {
  TIMED_SCOPE(t, "solveEquations");
  if (firstRun) {
    initializePARDISO(matrix);
  }
  firstRun = false;

  CHECK(numberOfEquations == matrix->getNumberOfRows());
  
	int phase;

  // initialize error code
	int error = 0; 

	// phase 22 is the numerical factorization
	phase = 22;
  int n = static_cast<int> (numberOfEquations);

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, matrix->getValuesArray(), 
      (int*) matrix->getIofeirArray(),
      (int*) matrix->getColumnsArray(),
			NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
  CHECK(error == 0) << "ERROR during numerical factorization. Error code = " << error;

	//Back substitution and iterative refinement
	phase = 33;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,	&n, matrix->getValuesArray(),
      (int*) matrix->getIofeirArray(),
      (int*) matrix->getColumnsArray(),
			NULL, &nrhs, iparm, &msglvl, rhs, unknowns, &error);

	CHECK(error == 0) << "ERROR during solution. Error code = " << error;
}

void PARDISO_equationSolver::initializePARDISO (math::SparseSymmetricMatrix* matrix) {
	for (uint16 i = 0; i < 64; i++) {
    iparm[i]=0;
  }

  for (uint16 i = 0; i < 64; i++) {
    pt[i]=0;
  }

	iparm[0] = 1; //no solver default
	iparm[1] = 2; //fill-in reordering from meris
	iparm[2] = MKL_Get_Max_Threads();
	iparm[3] = 0; //no iterative-direct algorithm
	iparm[4] = 0; //no user fill-in reducing permutation
	iparm[5] = 0; //write solution into x
	iparm[6] = 16; //default logical fortran unit number for output
	iparm[7] = 2; //max numbers of iterative refinement steps
	iparm[9] = 13; //pertrub the pivor elements with 1e-13
	iparm[10] = 1; //use nonsymmetric permutation  and scaling MPS
	iparm[13]=0; //output: number of perturbed pivots
	iparm[17]=-1; //output: number of nonzeros in the factor LU
	iparm[18]=-1; //output: MFLOPS for LU factorization
	iparm[19] = 0; //output: number of CG Iterations

  LOG_IF(!isSymmetric, FATAL) << "For now PARDISO_equationSolver doesn't support non-symmetric matrices";
	if (isPositive) {
    mtype = 2;
    LOG(INFO) << "EquationSolver will use positive symmetric solver";
  } else {
    LOG(INFO) << "EquationSolver will use non-positive symmetric solver";
    mtype = -2;
  }

  numberOfEquations = matrix->getNumberOfRows();
  int n = static_cast<int> (numberOfEquations);

  // initialize error code
	int error = 0; 

  int phase = 11;

  PARDISO(pt, &maxfct, &mnum, &mtype,&phase, &n, matrix->getValuesArray(),
     (int*) matrix->getIofeirArray(), 
     (int*) matrix->getColumnsArray(), 
			NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
  CHECK(error == 0) << "ERROR during symbolic factorization. Error code = " << error;
  LOG(INFO) << "Number of nonzeros in factors = " << iparm[17] << ", number of factorization MFLOPS = " << iparm[18];
}

void PARDISO_equationSolver::releasePARDISO () {
  int phase = -1;
  int n = static_cast<int> (numberOfEquations);

  // initialize error code
	int error = 0; 
	//Termination and release of memory
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,	&n, NULL, NULL, NULL, NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
  LOG_IF (error != 0, WARNING) << "ERROR during PARDISO termination. Error code = " << error;
}

} //namespace math

} //namespace nla3d

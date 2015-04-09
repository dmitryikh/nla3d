// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include <mkl.h>
#include "sys.h"

namespace nla3d {

// forward declarations
class FEStorage;
class EquationSolver;
class PostProcessor;

#define SOL_WAIT 1
#define SOL_SOLVING 2

#define SAE_CHOLESKY 1
#define SAE_BUNCH 2
#define SAE_DSS 3

/* класс содержит все механизмы для составления СЛАУ и Решения */
// How to refactor Solution class:
// 1.name FESolver
// 2. insert post processors here
// 3. be Transient solver ready (incorporate time, substeps, loadsteps - like in Ansys)
// 4. be ready for adaptive time stepping (be ready for restarts)
// 5. Able to solve linear problems, steady non-linear, transient, and buckling? / modal?.
//

class FESolver {
public:
  FESolver ();
  // attach FESolver to FEStorage instance
	void attachFEStorage (FEStorage *st);
	virtual void solve () = 0;
protected:
  FEStorage* storage = nullptr;
  EquationSolver* eqSolver = nullptr;
};

// particular realisation of FESolver for linear tasks
class LinearFESolver : public FESovler {
public:
  LinearFESolver ();
  void solve ();
};

// Iterative solver for Newton-Raphson procedure
// The convergence is controlled by mean deltaSolution (see solve procedure)
class NonlinearFESolver : public FESolver {
public:
  uint16 numberOfIterations = 20;
  uint16 numberOfLoadsteps = 10;

  double startTime = 0.0;
  double endTime = 1.0;

  double convergenceCriteria = 1.0e-3;

  void solve ();
protected:
};



class FESolver {
public:
	FESolver () {
		storage = NULL;
		status = SOL_WAIT;
		setInit();
		curCriteria = 0.0f;
		curLoadstep = 0;
		curIterat = 0;
		solver_type = SAE_CHOLESKY;
		if_solver_first_time = true;
	};

  // attach solver to FEStorage instance
	void attach (FEStorage *st);

	void run ();

	uint16 getqIterat ();
	uint16 getqLoadstep ();
	double getCriteria ();
	uint16 getStatus ();
	bool setqIterat (uint16 iter);
	bool setqLoadstep (uint16 ls);
	bool setCriteria (double cr);
	bool setInit ();
	
private:

	void *pt[64]; //Internal solver memory pointer pt
	bool if_solver_first_time;
	uint16 qLoadstep;	// кол-во шагов нагружений
	uint16 qIterat;		// кол-во итераций
	double Criteria;	// кинематический критерий сходимости
	uint16 status;
	uint16 solver_type;
	bool stopit;

	double curCriteria; // текущее значение критерия на данном шаге решения
	uint16 curLoadstep;
	uint16 curIterat;

	uint16 solve_wrap ();
	uint16 solveSAE_Bunch ();
	uint16 solveSAE_Cholesky ();
	uint16 solveSAE_DSS ();
	FEStorage* storage;
  // main process of solution
	void main_process ();
};

inline uint16 FESolver::getqIterat ()
{
	return qIterat;
}
inline uint16 FESolver::getqLoadstep ()
{
	return qLoadstep;
}
inline double FESolver::getCriteria ()
{
	return Criteria;
}
inline uint16 FESolver::getStatus ()
{
	return status;
}

} // namespace nla3d

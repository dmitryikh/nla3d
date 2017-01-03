// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "math/EquationSolver.h"

#include <Eigen/Dense>

namespace nla3d {

// forward declarations
class FEStorage;
class PostProcessor;

/* класс содержит все механизмы для составления СЛАУ и Решения */
// How to refactor Solution class:
// +1.name FESolver
// 2. insert post processors here
// +3. be Transient solver ready (incorporate time, substeps, loadsteps - like in Ansys)
// +4. be ready for adaptive time stepping (be ready for restarts)
// 5. be ready for restarts

class TimeControl {
public:
  uint16 getCurrentStep ();
  uint16 getNumberOfConvergedSteps ();

  uint16 getCurrentEquilibriumStep ();
  uint16 getTotalNumberOfEquilibriumSteps ();

  bool nextStep (double delta);
  void nextEquilibriumStep ();

  double getCurrentTime ();
  double getCurrentNormalizedTime ();
  double getEndTime ();
  double getStartTime ();
  double getCurrentTimeDelta ();
  double getCurrentNormalizedTimeDelta ();

  void setEndTime (double _endTime);
  void setStartTime (double _startTime);
protected:
  std::list<double> convergedTimeInstances;
  std::list<uint16> equilibriumSteps;
  uint16 currentEquilibriumStep = 0;
  uint16 totalNumberOfEquilibriumSteps = 0;
  double currentTimeDelta = 0.0;
  double endTime = 1.0;
  double startTime = 0.0;
  double currentTime = 0.0;
};

class FESolver {
public:
  FESolver ();
  ~FESolver ();
  // attach FESolver to FEStorage instance
	void attachFEStorage (FEStorage *st);
  void attachEquationSolver (math::EquationSolver *eq);
	virtual void solve () = 0;

	size_t getNumberOfPostProcessors ();
  // get an instance of PostProcessor by its number
  // _np >= 0
	PostProcessor& getPostProcessor(size_t _np);
  // The function stores PostProcessor in FESolver. Then FESolver takes care about 
  // dynamically allocated object - FESolver will free the memory.
	uint16 addPostProcessor (PostProcessor *pp);
  void deletePostProcessors ();
protected:
  FEStorage* storage = nullptr;
  math::EquationSolver* eqSolver = nullptr;

  // Array of PostProcessors. Here is a conception of PostProcessors that can do useful work every iteration
  // of the solution. For example here is a VtkProcessor class to write *.vtk file with model and solution data
  // (displacements, stresses and others). As Mpc, an instance of PostProcessor class is created outside
  // of FESolver. But with function addPostProcessor(..) FESolver take control on the PostProcessor instance.
  // They are deleted in ~FESolver() destructor.
	std::vector<PostProcessor*> postProcessors;
};

// particular realisation of FESolver for linear tasks
class LinearFESolver : public FESolver {
public:
  LinearFESolver();
  virtual void solve();
};

// Iterative solver for Newton-Raphson procedure
// The convergence is controlled by mean deltaSolution (see solve procedure)
class NonlinearFESolver : public FESolver {
public:
  NonlinearFESolver();

  TimeControl timeControl;
  uint16 numberOfIterations = 20;
  uint16 numberOfLoadsteps = 10;

  double convergenceCriteria = 1.0e-3;

  virtual void solve();
protected:
  double calculateCriteria();
};


// use Newmark scheme
class LinearTransientFESolver : public FESolver {
public:
  LinearTransientFESolver();

  TimeControl timeControl;
  uint16 numberOfIterations = 1;
  uint16 numberOfTimesteps = 100;

  double alpha = 0.25; // trapezoidal rule by default
  double delta = 0.5;

  double time0 = 0.0;
  double time1 = 1.0;
  double initValue = 0.0;

  virtual void solve ();

  double a0, a1, a2, a3, a4, a5, a6, a7;
};

} // namespace nla3d

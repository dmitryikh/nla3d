// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "math/Vec.h"
#include "math/EquationSolver.h"
#include "FEStorage.h"
#include "PostProcessor.h"

#include <Eigen/Dense>

namespace nla3d {

using namespace math;

class FEStorage;
class PostProcessor;

// TimeControl class helps to control time stepping and equilibrium iterations.
// Before solution one can setup setStartTime() and setEndTime() and then perform time stepping with
// nextStep(delta) where delta is solver specific time step. For every time step many equilibrium
// iterations could be performed by nextEquilibriumStep(). All intermediate time steps are stored in
// convergedTimeInstances, for every time steps number of equilibrium steps are stored in
// equilibriumSteps.
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

// The abstract class that represents the FE solver. This class use information and methods from
// FEStorage in order to apply particular solution scheme. For example, here are already some of
// them: solver for linear steady FE model with concentrated loads and fixed DoFs (LinearFESolver),
// solver for non-linear quasi-steady FE model based on incremental approach (NonlinearFESolver),
// solver for linear transient FE model based on Newmark time-integration scheme
// (LinearTransientFESolver).
//
// FESolver keeps a pointer to FEStorage to receive assembled global equation system, and a pointer
// to math::EquationSolver to solve linear system of equations.
//
// FESovler manages PostProcessors objects in order to call them on every time step to do some
// useful routines (write res files, log reaction forces, etc.)
//
// FESolver by default contains arrays of loadBC and fixBC structures in order to apply constant BC
// on DoFs (concentrated load and fixed DoF value). But particular realizations of FESolver can
// incorporate time-depended BC along with or instead of loadBC and fixBC.
//
// FESolver by mean of initSolutionData() keeps pointers to FEStorage matrices and vectors. matK for
// stiffness matrix, matC for dumping matrix, matM for inertia matrix. vecU - for DoF values vector
// with references to particular parts of it: vecUc(part for constrained DoFs), vecUs(part for
// unknown DoFs), vecUl(part for Lagrangian lambdas from Mpc equations), vecUsl(combination of vecUs and
// vecUl). The same for first and second time derivatives: vecDU, vecDDU.
// Also there are references to RHS parts: vecF - for FE internal loads. vecFl has special treatment
// as RHS of Mpc equations. vecR - for external concentrated loads. vecRc is treated as reaction
// forces of fixation for constrained DoFs, vecRs is treated as external concetrated loads and vecRl
// should be zero all times. All references above are just easy way to access FEStorage internal
// data. Actually all this vectors and matrices ale allocated inside of FEStorage. By performing
// FEStorage::assembleGlobalEqMatrices() matK/C/M, vecF are assembled with particular numbers based
// on FE element formulation and on Mpc equations. Other entities (vecU, vecDU, vecDDU, vecR) are
// fully under FESovler control. FESolver is responsible on timely updating this values.
//
class FESolver {
  public:
    FESolver ();
    virtual ~FESolver ();

    void attachFEStorage(FEStorage *st);
    void attachEquationSolver(math::EquationSolver *eq);

    // main method were all solution scheme specific routines are performed
    virtual void solve() = 0;

    size_t getNumberOfPostProcessors();
    // get an instance of PostProcessor by its number
    // _np >= 0
    PostProcessor& getPostProcessor(size_t _np);
    // The function stores PostProcessor pointer in FESolver. FESovler will free the memory by itself.
    uint16 addPostProcessor(PostProcessor *pp);
    void deletePostProcessors();

    // run FEStorage procedures for initialization solution data structures and map matK, matC, matM,
    // vecU, vecDU, vecDDU, vecR, vecF pointer on FEStorage's allocated structures.
    void initSolutionData();

    // The function applies boundary conditions stored in `loads` and `fixs` arrays.
    // `time` should be normalized(in range [0.0; 1.0])
    void applyBoundaryConditions(double time);

    // constrained DoFs has special treatment in FEStorage, that's why before initialization of
    // solution data we need to tell FEStorage which DoFs is constrained. setConstrainedDofs() does
    // this work based on `fixs` array.
    void setConstrainedDofs();

    // add DoF fixation (constraint) boundary condition 
    void addFix(int32 n, Dof::dofType dof, const double value = 0.0);
    // add DoF load (force) boundary condition 
    void addLoad(int32 n, Dof::dofType dof, const double value = 0.0);
    
    // for debug purpose:
    // dump matrices matK, matC, matM and vectors vecF, vecR
    void dumpMatricesAndVectors(std::string filename);
    // read matrices matK, matC, matM and vectors vecF, vecR from file and compare them with
    // solution instances
    void compareMatricesAndVectors(std::string filename, double th = 1.0e-9);
  protected:
    FEStorage* storage = nullptr;
    math::EquationSolver* eqSolver = nullptr;

    // Array of PostProcessors. Here is a conception of PostProcessors that can do useful work every
    // iteration of the solution. For example here is a VtkProcessor class to write *.vtk file with
    // model and solution data (displacements, stresses and others). An instance of PostProcessor
    // class is created outside of FESolver. But with function addPostProcessor(..) FESolver take
    // control on the PostProcessor instance.  The memory released in deletePostProcessors().
    std::vector<PostProcessor*> postProcessors;

    // the references on FEStorage FE data structures which is frequently used by FESolver. This
    // references make it easy to operate with important FE data structures. 
    BlockSparseSymMatrix<2>* matK = nullptr;
    BlockSparseSymMatrix<2>* matC = nullptr;
    BlockSparseSymMatrix<2>* matM = nullptr;

    // vecU is a reference on a whole DoF values vector, but vecUc, vecUs, vecUl, vecUsl are
    // references on parts of the full vector. Some times it's handy to operate with the particular
    // part only.
    dVec vecU;
    dVec vecUc;
    dVec vecUs;
    dVec vecUl;
    dVec vecUsl;

    // for first time derivatives
    dVec vecDU;
    dVec vecDUc;
    dVec vecDUs;
    dVec vecDUl;
    dVec vecDUsl;

    // for second time derivatives
    dVec vecDDU;
    dVec vecDDUc;
    dVec vecDDUs;
    dVec vecDDUl;
    dVec vecDDUsl;

    // for internal element loads
    dVec vecF;
    dVec vecFc;
    dVec vecFs;
    dVec vecFl;
    dVec vecFsl;

    // for external loads
    dVec vecR;
    dVec vecRc;
    dVec vecRs;
    dVec vecRl;
    dVec vecRsl;

    // List of concentrated load boundary conditions (addition to RHS of global eq. system)
    std::list<loadBC> loads;

    // List of prescribed values for DoFs.
    // NOTE: DoFs with fixed values are treated in nla3d in a special manner: this DoFs are eliminated from
    // global solve system and then reaction loads (vecRc vector) of a such DoFs are found.
    std::list<fixBC> fixs;
};


// particular realization of FESolver for linear tasks
class LinearFESolver : public FESolver {
  public:
    LinearFESolver();
    virtual void solve();
};

// Iterative solver for Full Newton-Raphson procedure
// The convergence is controlled by mean increment of DoF values
class NonlinearFESolver : public FESolver {
  public:
    NonlinearFESolver();

    TimeControl timeControl;
    uint16 numberOfIterations = 20;
    uint16 numberOfLoadsteps = 10;

    double convergenceCriteria = 1.0e-3;

    virtual void solve();
  protected:
    double calculateCriteria(dVec& delta);
};


// Solver for time integration of linear systems: M * DDU + C * DU + K * U = F + R
// use Newmark scheme (based on section 9.2.4. Bathe K.J., Finite Element Procedures, 1997)
class LinearTransientFESolver : public FESolver {
  public:
    LinearTransientFESolver();

    uint16 numberOfTimesteps = 100;

    double alpha = 0.25; // trapezoidal rule by default
    double delta = 0.5;

    // start time
    double time0 = 0.0;
    // end time
    double time1 = 1.0;
    // initial values for vecUc
    double initValue = 0.0;

    virtual void solve ();

    double a0, a1, a2, a3, a4, a5, a6, a7;
};

} // namespace nla3d

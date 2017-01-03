// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "FESolver.h"

#include "FEStorage.h"
#include "PostProcessor.h"

namespace nla3d {

using namespace ::nla3d::math;


uint16 TimeControl::getCurrentStep () {
  return static_cast<uint16> (convergedTimeInstances.size()+1);
}

uint16 TimeControl::getNumberOfConvergedSteps () {
  return static_cast<uint16> (convergedTimeInstances.size());
}

uint16 TimeControl::getCurrentEquilibriumStep () {
  return currentEquilibriumStep;
}

uint16 TimeControl::getTotalNumberOfEquilibriumSteps () {
  return totalNumberOfEquilibriumSteps; 
}

bool TimeControl::nextStep (double delta) {
  if (currentEquilibriumStep > 0) {
    equilibriumSteps.push_back (currentEquilibriumStep);
    convergedTimeInstances.push_back (currentTime);
  }
  if (currentTime >= endTime) {
    return false;    
  }
  // if delta bigger than endTime-currentTime
  if (currentTime + delta >= endTime) {
    currentTimeDelta = endTime - currentTime;
  } else if (currentTime + 2*delta >= endTime) {
    currentTimeDelta = (endTime - currentTime)/2.0;
  } else {
    currentTimeDelta = delta;
  }
  DCHECK (currentTimeDelta > 0.0);
  currentTime += currentTimeDelta;
  currentEquilibriumStep = 0;
  LOG(INFO) << "***** Loadstep = " << getCurrentStep() << ", Time = "
      << getCurrentTime() << " of " << getEndTime();
  return true;
}

void TimeControl::nextEquilibriumStep () {
  currentEquilibriumStep++;
  totalNumberOfEquilibriumSteps++; 
  LOG(INFO) << "***** Equilibrium iteration = " << currentEquilibriumStep
      << ", Cumulutive iterations = " << totalNumberOfEquilibriumSteps;
}

double TimeControl::getCurrentTime () {
  return currentTime;
}

double TimeControl::getCurrentNormalizedTime () {
  return (currentTime - startTime) / (endTime - startTime);
}

double TimeControl::getEndTime () {
  return endTime;
}

double TimeControl::getStartTime () {
  return startTime;
}

double TimeControl::getCurrentTimeDelta () {
  if (currentEquilibriumStep == 1) {
    return currentTimeDelta;
  } else {
    return 0.0;
  }
}

double TimeControl::getCurrentNormalizedTimeDelta () {
  if (currentEquilibriumStep == 1) {
    return currentTimeDelta / (endTime - startTime);
  } else {
    return 0.0;
  }
}

void TimeControl::setEndTime (double _endTime) {
  LOG_IF(currentTime > 0.0, ERROR) << "Trying to set end time = " << _endTime
      << " when solution is running (current time = " << currentTime << ")";
  endTime = _endTime;
}

void TimeControl::setStartTime (double _startTime) {
  LOG_IF(currentTime > 0.0, ERROR) << "Trying to set start time = " << _startTime
      << " when solution is running (current time = " << currentTime << ")";
  startTime = _startTime;
}


FESolver::FESolver () {
  eqSolver = defaultEquationSolver;
}

FESolver::~FESolver () {
  deletePostProcessors();
}

void FESolver::attachFEStorage (FEStorage *st) {
	if (storage) {
    LOG(WARNING) << "FEStorage already is attached. The old one will be dropped.";
  }
	storage = CHECK_NOTNULL(st);
}

size_t FESolver::getNumberOfPostProcessors () {
  return postProcessors.size();
}

// numbering from 0
PostProcessor& FESolver::getPostProcessor(size_t _np) {
	CHECK (_np < getNumberOfPostProcessors());
	return *postProcessors[_np];
}

uint16 FESolver::addPostProcessor (PostProcessor *pp) {
	CHECK_NOTNULL (pp);
	uint16 num = static_cast<uint16> (this->postProcessors.size()+1);
	pp->nPost_proc = num;
	postProcessors.push_back(pp);
	return num;
}

void FESolver::deletePostProcessors () {
	for (size_t i = 0; i < postProcessors.size(); i++) {
		delete postProcessors[i];
  }
  postProcessors.clear();
}

void FESolver::attachEquationSolver (math::EquationSolver *eq) {
	if (eqSolver) {
    LOG(WARNING) << "EquationSolver already is attached. The old one will be dropped.";
  }
	eqSolver = CHECK_NOTNULL(eq);
}


LinearFESolver::LinearFESolver () : FESolver () {

}

void LinearFESolver::solve () {
  TIMED_SCOPE(timer, "solution");
	LOG(INFO) << "Start the solution process";
  CHECK_NOTNULL(storage);
  CHECK_NOTNULL(eqSolver);

  // setup matrix properties for EquationSolver 
  eqSolver->setSymmetric(true);
  eqSolver->setPositive(false);

  bool success = storage->initializeSolutionData();
  LOG_IF(success != true, FATAL) << "initializeSolutionData is failed";

	for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
		postProcessors[i]->pre();
  }

  storage->assembleGlobalEqMatrices();
  storage->applyBoundaryConditions(1.0, 1.0);
  // solve equation system
  eqSolver->solveEquations (storage->getK(), storage->getF(), storage->getGlobalEqUnknowns());

  storage->updateSolutionResults();

  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->process (1);
  }

  LOG(INFO) << "***** SOLVED *****";

  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->post (1);
  }
}


NonlinearFESolver::NonlinearFESolver() : FESolver () {

}

void NonlinearFESolver::solve() {
  TIMED_SCOPE(timer, "solution");
	LOG(INFO) << "Start the solution process";
  CHECK_NOTNULL(storage);
  CHECK_NOTNULL(eqSolver);

  // setup matrix properties for EquationSolver 
  eqSolver->setSymmetric(true);
  eqSolver->setPositive(false);

  bool success = storage->initializeSolutionData();
  LOG_IF(success != true, FATAL) << "initializeSolutionData is failed";

	for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
		postProcessors[i]->pre();
  }

	double currentCriteria = 0.0;
  double timeDelta = (timeControl.getEndTime() - timeControl.getStartTime()) / numberOfLoadsteps;

  while (timeControl.nextStep(timeDelta)) {
    bool converged = false;
    for (;;) {
      timeControl.nextEquilibriumStep();

      storage->assembleGlobalEqMatrices();
			storage->applyBoundaryConditions(timeControl.getCurrentNormalizedTime(),
          timeControl.getCurrentNormalizedTimeDelta());

			// solve equation system
      eqSolver->solveEquations (storage->getK(), storage->getF(), storage->getGlobalEqUnknowns());

      // for now error in equation solver leads to FATAL error. In this case we don't need success
      // check on this level..
      //LOG_IF(success != true, FATAL) << "equation solver failed";

			storage->updateSolutionResults();

      // calculate convergence criteria
      currentCriteria = calculateCriteria();

			if (currentCriteria < convergenceCriteria) {
        converged = true;
				break;
			}
      LOG_IF(currentCriteria > 1.0e6 || std::isnan(currentCriteria), FATAL) << "The solution is diverged!";
		}//iterations

    LOG_IF(!converged, FATAL) << "The solution is not converged with"
        << timeControl.getCurrentEquilibriumStep() << " equilibrium iterations";
		LOG(INFO) << "Loadstep " << timeControl.getCurrentStep() << " completed with " << timeControl.getCurrentEquilibriumStep();

// TODO: figure out why TIMED_BLOCK doesn't work here..
//    TIMED_BLOCK(t, "PostProcessor::process") {
        for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
          postProcessors[i]->process (timeControl.getCurrentStep());
        }
//    }
	} //loadsteps
  LOG(INFO) << "***** SOLVED *****";

//  TIMED_BLOCK(t, "PostProcessor::post") {
      for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
        postProcessors[i]->post (timeControl.getCurrentStep());
      }
//  }
}

double NonlinearFESolver::calculateCriteria() {
  double curCriteria = 0.0;
  for (uint32 i=0; i < storage->getNumberOfUnknownDofs(); i++) {
    curCriteria += fabs(storage->getGlobalEqUnknowns()[i]);
  }
  curCriteria /= storage->getNumberOfUnknownDofs();

  curCriteria = curCriteria / convergenceCriteria;
  LOG(INFO) << "Solution delta criteria = " << curCriteria;
  return curCriteria;
}


LinearTransientFESolver::LinearTransientFESolver() : FESolver () {

}

void LinearTransientFESolver::solve() {
  TIMED_SCOPE(timer, "solution");
	LOG(INFO) << "Start the solution process";
  CHECK_NOTNULL(storage);
  CHECK_NOTNULL(eqSolver);

  // setup matrix properties for EquationSolver 
  eqSolver->setSymmetric(true);
  eqSolver->setPositive(false);

  storage->setTransient(true);

  double dt = (time1 - time0) / numberOfTimesteps;
  double curTime = time0 + dt;
  uint16 curTimestep = 1;
  // initialize parameters of Newmark procedure
  a0 = 1.0 / (alpha * dt * dt);
  a1 = delta / (alpha * dt);
  a2 = 1.0 / (alpha * dt);
  a3 = 1.0 / (2.0 * alpha) - 1.0;
  a4 = delta / alpha - 1.0;
  a5 = dt / 2.0 * (delta / alpha - 2.0);
  a6 = dt * (1.0 - delta);
  a7 = delta * dt;

  bool success = storage->initializeSolutionData();
  LOG_IF(success != true, FATAL) << "initializeSolutionData is failed";

	for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
		postProcessors[i]->pre();
  }

  storage->assembleGlobalEqMatrices();
  storage->applyBoundaryConditions(1.0, 1.0);

  auto matK = storage->getK();
  Eigen::Map<Eigen::VectorXd> eK(matK->getValuesArray(), matK->nValues());

  auto matM = storage->getM();
  Eigen::Map<Eigen::VectorXd> eM(matM->getValuesArray(), matM->nValues());

  auto matC = storage->getC();
  Eigen::Map<Eigen::VectorXd> eC(matC->getValuesArray(), matC->nValues());

  uint32 nEq = matK->nRows();

  Eigen::Map<Eigen::VectorXd> vecF(storage->getF(), nEq);

  math::SparseSymMatrix matKmod(matK->getSparsityInfo());
  Eigen::Map<Eigen::VectorXd> eKmod(matKmod.getValuesArray(), matKmod.nValues());

  // thanks to Eigen map to my values memory
  eKmod = a0 * eM + a1 * eC + eK;

  // factorize eKmode just once as far sa we have a deal with linear system
  eqSolver->factorizeEquations(&matKmod);

  Eigen::VectorXd vecFmod = Eigen::VectorXd::Zero(nEq);

  // NOTE: initial conditions are zeros!
  Eigen::VectorXd U = Eigen::VectorXd::Zero(nEq);
  for (auto i = 0; i < nEq; i++) 
    U[i] = initValue;
  Eigen::VectorXd Unext = Eigen::VectorXd::Zero(nEq);
  Eigen::VectorXd Udot = Eigen::VectorXd::Zero(nEq);
  Eigen::VectorXd Udotnext = Eigen::VectorXd::Zero(nEq);
  Eigen::VectorXd Udotdot = Eigen::VectorXd::Zero(nEq);
  Eigen::VectorXd Udotdotnext = Eigen::VectorXd::Zero(nEq);

  Eigen::VectorXd tmp1 = Eigen::VectorXd::Zero(nEq);
  Eigen::VectorXd tmp2 = Eigen::VectorXd::Zero(nEq);

  while (curTime <= time1) {
    for (;;) {
      // restore rhs into vecFmod
      tmp1 = a0 * U + a2 * Udot + a3 * Udotdot;
      tmp2 = a1 * U + a4 * Udot + a5 * Udotdot;
      for (uint32 i = 0; i < nEq; i++)
        vecFmod[i] = vecF[i] + matM->mult_vec_i(tmp1.data(), i+1) + matC->mult_vec_i(tmp2.data(), i + 1);

			// solve equation system
      eqSolver->substituteEquations(&matKmod, &vecFmod[0], &Unext[0]);

      break;
		}//iterations

    // one more Kludge 
    double* stU = storage->getGlobalEqUnknowns();
    for (uint32 i = 0; i < nEq; i++) {
      *stU = Unext[i];
      stU++;
    }

    storage->updateTimestepResults();

    // restore derivatives
    Udotdotnext = a0 * (Unext - U) - a2 * Udot - a3 * Udotdot;
    Udotnext = Udot + a6 * Udotdotnext + a7 * Udotdot;

    U = Unext;
    Udot = Udotnext;
    Udotdot = Udotdotnext;

		LOG(INFO) << "Time " << curTime << " completed";

    for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
      postProcessors[i]->process(curTimestep);
    }

    curTimestep++;
    curTime += dt;
	} //timesteps
  LOG(INFO) << "***** SOLVED *****";

  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->post (timeControl.getCurrentStep());
  }
}

} // namespace nla3d

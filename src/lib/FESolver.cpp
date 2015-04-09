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

  storage->assembleGlobalEqMatrix();
  storage->applyBoundaryConditions(1.0, 1.0);
  // solve equation system
  eqSolver->solveEquations (storage->getGlobalEqMatrix(), storage->getGlobalEqRhs(), storage->getGlobalEqUnknowns());

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

      storage->assembleGlobalEqMatrix();
			storage->applyBoundaryConditions(timeControl.getCurrentNormalizedTime(),
          timeControl.getCurrentNormalizedTimeDelta());

			// solve equation system
      eqSolver->solveEquations (storage->getGlobalEqMatrix(), storage->getGlobalEqRhs(), storage->getGlobalEqUnknowns());

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


} // namespace nla3d

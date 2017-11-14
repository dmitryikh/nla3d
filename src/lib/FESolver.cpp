// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "FESolver.h"

namespace nla3d {

using namespace ::nla3d::math;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// TimeControl
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
uint16 TimeControl::getCurrentStep() {
  return static_cast<uint16> (convergedTimeInstances.size()+1);
}


uint16 TimeControl::getNumberOfConvergedSteps() {
  return static_cast<uint16> (convergedTimeInstances.size());
}


uint16 TimeControl::getCurrentEquilibriumStep() {
  return currentEquilibriumStep;
}


uint16 TimeControl::getTotalNumberOfEquilibriumSteps() {
  return totalNumberOfEquilibriumSteps; 
}


bool TimeControl::nextStep(double delta) {
  if (currentEquilibriumStep > 0) {
    equilibriumSteps.push_back(currentEquilibriumStep);
    convergedTimeInstances.push_back(currentTime);
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


void TimeControl::nextEquilibriumStep() {
  currentEquilibriumStep++;
  totalNumberOfEquilibriumSteps++; 
  LOG(INFO) << "***** Equilibrium iteration = " << currentEquilibriumStep
      << ", Cumulutive iterations = " << totalNumberOfEquilibriumSteps;
}


double TimeControl::getCurrentTime() {
  return currentTime;
}


double TimeControl::getCurrentNormalizedTime() {
  return (currentTime - startTime) / (endTime - startTime);
}


double TimeControl::getEndTime() {
  return endTime;
}


double TimeControl::getStartTime() {
  return startTime;
}


double TimeControl::getCurrentTimeDelta() {
  if (currentEquilibriumStep == 1) {
    return currentTimeDelta;
  } else {
    return 0.0;
  }
}


double TimeControl::getCurrentNormalizedTimeDelta() {
  if (currentEquilibriumStep == 1) {
    return currentTimeDelta / (endTime - startTime);
  } else {
    return 0.0;
  }
}


void TimeControl::setEndTime(double _endTime) {
  LOG_IF(currentTime > 0.0, ERROR) << "Trying to set end time = " << _endTime
      << " when solution is running (current time = " << currentTime << ")";
  endTime = _endTime;
}


void TimeControl::setStartTime(double _startTime) {
  LOG_IF(currentTime > 0.0, ERROR) << "Trying to set start time = " << _startTime
      << " when solution is running (current time = " << currentTime << ")";
  startTime = _startTime;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// FESolver
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
FESolver::FESolver() {
  eqSolver = defaultEquationSolver;
}


FESolver::~FESolver() {
  deletePostProcessors();
}


void FESolver::attachFEStorage(FEStorage *st) {
  if (storage) {
    LOG(WARNING) << "FEStorage already is attached. The old one will be dropped.";
  }
  storage = CHECK_NOTNULL(st);
}


size_t FESolver::getNumberOfPostProcessors() {
  return postProcessors.size();
}

// numbering from 0
PostProcessor& FESolver::getPostProcessor(size_t _np) {
  CHECK (_np < getNumberOfPostProcessors());
  return *postProcessors[_np];
}


uint16 FESolver::addPostProcessor(PostProcessor *pp) {
  CHECK_NOTNULL (pp);
  uint16 num = static_cast<uint16>(this->postProcessors.size()+1);
  pp->nPost_proc = num;
  postProcessors.push_back(pp);
  return num;
}


void FESolver::deletePostProcessors() {
  for (size_t i = 0; i < postProcessors.size(); i++) {
    delete postProcessors[i];
  }
  postProcessors.clear();
}


void FESolver::attachEquationSolver(math::EquationSolver *eq) {
  if (eqSolver) {
    LOG(WARNING) << "EquationSolver already is attached. The old one will be dropped.";
  }
  eqSolver = CHECK_NOTNULL(eq);
}

void FESolver::initSolutionData() {
  storage->initSolutionData();

  // make easy access to global system of equations entities
  // pointer to stiff. matrix
  matK = storage->getK();

  // pointers to DoF values vector and it parts (c - constrained DoFs, s - to be solved (unknown)
  // DoFs, l - mpc's lambdas
  vecU.reinit(*(storage->getU()), 0, storage->nConstrainedDofs() + storage->nUnknownDofs() + storage->nMpc());
  vecUc.reinit(*(storage->getU()), 0, storage->nConstrainedDofs());
  vecUs.reinit(*(storage->getU()), storage->nConstrainedDofs(), storage->nUnknownDofs());
  vecUl.reinit(*(storage->getU()), storage->nConstrainedDofs() + storage->nUnknownDofs(), storage->nMpc());
  vecUsl.reinit(*(storage->getU()), storage->nConstrainedDofs(), storage->nUnknownDofs() + storage->nMpc());

  vecF.reinit(*(storage->getF()), 0, storage->nConstrainedDofs() + storage->nUnknownDofs() + storage->nMpc());
  vecFc.reinit(*(storage->getF()), 0, storage->nConstrainedDofs());
  vecFs.reinit(*(storage->getF()), storage->nConstrainedDofs(), storage->nUnknownDofs());
  vecFl.reinit(*(storage->getF()), storage->nConstrainedDofs() + storage->nUnknownDofs(), storage->nMpc());
  vecFsl.reinit(*(storage->getF()), storage->nConstrainedDofs(), storage->nUnknownDofs() + storage->nMpc());

  vecR.reinit(*(storage->getR()), 0, storage->nConstrainedDofs() + storage->nUnknownDofs() + storage->nMpc());
  vecRc.reinit(*(storage->getR()), 0, storage->nConstrainedDofs());
  vecRs.reinit(*(storage->getR()), storage->nConstrainedDofs(), storage->nUnknownDofs());
  vecRl.reinit(*(storage->getR()), storage->nConstrainedDofs() + storage->nUnknownDofs(), storage->nMpc());
  vecRsl.reinit(*(storage->getR()), storage->nConstrainedDofs(), storage->nUnknownDofs() + storage->nMpc());

  if (storage->isTransient()) {
    matC = storage->getC();
    matM = storage->getM();

    vecDU.reinit(*(storage->getDU()), 0, storage->nConstrainedDofs() + storage->nUnknownDofs() + storage->nMpc());
    vecDUc.reinit(*(storage->getDU()), 0, storage->nConstrainedDofs());
    vecDUs.reinit(*(storage->getDU()), storage->nConstrainedDofs(), storage->nUnknownDofs());
    vecDUl.reinit(*(storage->getDU()), storage->nConstrainedDofs() + storage->nUnknownDofs(), storage->nMpc());
    vecDUsl.reinit(*(storage->getDU()), storage->nConstrainedDofs(), storage->nUnknownDofs() + storage->nMpc());

    vecDDU.reinit(*(storage->getDDU()), 0, storage->nConstrainedDofs() + storage->nUnknownDofs() + storage->nMpc());
    vecDDUc.reinit(*(storage->getDDU()), 0, storage->nConstrainedDofs());
    vecDDUs.reinit(*(storage->getDDU()), storage->nConstrainedDofs(), storage->nUnknownDofs());
    vecDDUl.reinit(*(storage->getDDU()), storage->nConstrainedDofs() + storage->nUnknownDofs(), storage->nMpc());
    vecDDUsl.reinit(*(storage->getDDU()), storage->nConstrainedDofs(), storage->nUnknownDofs() + storage->nMpc());
  }
}


void FESolver::setConstrainedDofs() {
  for (auto& fix : fixs) {
    // TODO: now support only nodal dofs..
    storage->setConstrainedNodeDof(fix.node, fix.node_dof);
  }
}


void FESolver::applyBoundaryConditions(double time) {
  TIMED_SCOPE(t, "applyBoundaryConditions");
  LOG(INFO) << "Applying boundary conditions.. (" << loads.size() << " nodal loads and "
       << fixs.size() << "nodal fixations)";  

  // fill nodal loads
  for (auto& load : loads) {
    storage->addValueR(load.node, load.node_dof, load.value * time);
  }

  // fill nodal displacements (kinematic fixs)
  for (auto& fix : fixs) {
    //TODO: now support only nodal dofs..
    uint32 eq_num = storage->getNodeDofEqNumber(fix.node, fix.node_dof);
    // To be sure that constrained DoF lays in numberOfConstrainedDofs part
    assert (eq_num - 1 < storage->nConstrainedDofs());
    vecUc[eq_num - 1] = fix.value * time;
  }
}


void FESolver::addFix(int32 n, Dof::dofType dof, const double value) {
  // NOTE: we believes that all fixBC are distinct!
  //       Do not pass to it the same BC twice!
  fixs.push_back(fixBC(n, dof, value));
}


void FESolver::addLoad(int32 n, Dof::dofType dof, const double value) {
  loads.push_back(loadBC(n, dof, value));
}


void FESolver::dumpMatricesAndVectors(std::string filename) {
  std::ofstream out(filename);
  matK->block(1)->writeCoordinateTextFormat(out);
  matK->block(1, 2)->writeCoordinateTextFormat(out);
  matK->block(2)->writeCoordinateTextFormat(out);
  if (storage->isTransient()) {
    matC->block(1)->writeCoordinateTextFormat(out);
    matC->block(1, 2)->writeCoordinateTextFormat(out);
    matC->block(2)->writeCoordinateTextFormat(out);
    matM->block(1)->writeCoordinateTextFormat(out);
    matM->block(1, 2)->writeCoordinateTextFormat(out);
    matM->block(2)->writeCoordinateTextFormat(out);
  }

  vecF.writeTextFormat(out);
  vecR.writeTextFormat(out);

  out.close();
}


void FESolver::compareMatricesAndVectors(std::string filename, double th) {
  std::ifstream in(filename);
  SparseSymMatrix K1;
  K1.readCoordinateTextFormat(in);
  CHECK(matK->block(1)->compare(K1, th));
  SparseMatrix K12;
  K12.readCoordinateTextFormat(in);
  CHECK(matK->block(1, 2)->compare(K12, th));
  SparseSymMatrix K2;
  K2.readCoordinateTextFormat(in);
  CHECK(matK->block(2)->compare(K2, th));

  if (storage->isTransient()) {
    SparseSymMatrix C1;
    C1.readCoordinateTextFormat(in);
    CHECK(matC->block(1)->compare(C1, th));
    SparseMatrix C12;
    C12.readCoordinateTextFormat(in);
    CHECK(matC->block(1, 2)->compare(C12, th));
    SparseSymMatrix C2;
    C2.readCoordinateTextFormat(in);
    CHECK(matC->block(2)->compare(C2, th));

    SparseSymMatrix M1;
    M1.readCoordinateTextFormat(in);
    CHECK(matM->block(1)->compare(M1, th));
    SparseMatrix M12;
    M12.readCoordinateTextFormat(in);
    CHECK(matM->block(1, 2)->compare(M12, th));
    SparseSymMatrix M2;
    M2.readCoordinateTextFormat(in);
    CHECK(matM->block(2)->compare(M2, th));
  }

  dVec fvecF;
  fvecF.readTextFormat(in);
  CHECK(vecF.compare(fvecF, th));

  dVec fvecR;
  fvecR.readTextFormat(in);
  in.close();
  CHECK(vecR.compare(fvecR, th));

  in.close();
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// LinearFESolver
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
LinearFESolver::LinearFESolver() : FESolver() {

}

void LinearFESolver::solve () {
  TIMED_SCOPE(timer, "solution");
  LOG(INFO) << "Start the solution process";
  CHECK_NOTNULL(storage);
  CHECK_NOTNULL(eqSolver);

  // setup matrix properties for EquationSolver 
  eqSolver->setSymmetric(true);
  eqSolver->setPositive(false);

  // This is right procedures to init solution infrmation in FEStorage:
  // 1. Register all DoFs that will be used in solution
  storage->initDofs();
  // 2. tell FEStorage which DoFs will be fixed (constrained)
  setConstrainedDofs();
  // 3. Perform global system eq. numbering (first goes constrained DoFs, then unknown, Mpc
  // equations are last ones)
  storage->assignEquationNumbers();
  // 4. Allocate memory, initialize matrices, assign pointer on this data.
  initSolutionData();

  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->pre();
  }

  vecR.zero();

  storage->assembleGlobalEqMatrices();
  applyBoundaryConditions(1.0);

  // we need to calculate RHS for unknowns Dofs and Mps eq.
  dVec rhs(storage->nUnknownDofs() + storage->nMpc());

  rhs.zero();
  rhs += vecFsl;
  rhs += vecRsl;
  // need to take into account elimination of constrained DoFs:
  matBTVprod(*(matK->block(1,2)), vecUc, -1.0, rhs);

  // solve equation system
  eqSolver->solveEquations(matK->block(2), rhs.ptr(), vecUsl.ptr());

  // restore reaction loads for constrained DoFs.
  vecRc.zero();
  matBVprod(*(matK->block(1)), vecUc, 1.0, vecRc);
  matBVprod(*(matK->block(1,2)), vecUsl, 1.0, vecRc);
  vecRc -= vecFc;

  // update results for elements
  storage->updateResults();

  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->process (1);
  }

  LOG(INFO) << "***** SOLVED *****";

  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->post (1);
  }
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// NonlinearFESolver
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
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

  storage->initDofs();
  setConstrainedDofs();
  storage->assignEquationNumbers();
  initSolutionData();

  dVec rhs(storage->nUnknownDofs() + storage->nMpc());
  // need to store obtained constrained DoFs values obtained on previous equilibrium step in order
  // to compute deltaUc for incremental approach
  dVec Ucprev(storage->nConstrainedDofs());
  Ucprev.zero();
  // deltaU vector of unknowns for incremental approach
  dVec deltaU(storage->nConstrainedDofs() + storage->nUnknownDofs() + storage->nMpc());
  dVec deltaUc(deltaU, 0, storage->nConstrainedDofs());
  dVec deltaUsl(deltaU, storage->nConstrainedDofs(), storage->nUnknownDofs() + storage->nMpc());
  dVec deltaUs(deltaU, storage->nConstrainedDofs(), storage->nUnknownDofs());
  dVec deltaUl(deltaU, storage->nConstrainedDofs() + storage->nUnknownDofs(), storage->nMpc());

  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->pre();
  }

  double currentCriteria = 0.0;
  double timeDelta = (timeControl.getEndTime() - timeControl.getStartTime()) / numberOfLoadsteps;

  while (timeControl.nextStep(timeDelta)) {
    bool converged = false;
    for (;;) {
      timeControl.nextEquilibriumStep();
      vecR.zero();
      storage->assembleGlobalEqMatrices();
      applyBoundaryConditions(timeControl.getCurrentNormalizedTime());

      deltaUc = vecUc - Ucprev;
      rhs.zero();
      rhs += vecFsl;
      rhs += vecRsl;
      // nConstr x (nUnknown + nMPC)
      matBTVprod(*(matK->block(1,2)), deltaUc, -1.0, rhs);

      // solve equation system
      eqSolver->solveEquations(matK->block(2), rhs.ptr(), deltaUsl.ptr());

      // restore DoF values from increments
      vecUs += deltaUs;
      vecUl = deltaUl;

      // restore constrained DoFs reactions
      vecRc.zero();
      matBVprod(*(matK->block(1)), deltaUc, 1.0, vecRc);
      matBVprod(*(matK->block(1,2)), deltaUsl, 1.0, vecRc);
      vecRc -= vecFc;

      Ucprev = vecUc;

      storage->updateResults();

      // calculate convergence criteria
      currentCriteria = calculateCriteria(deltaUs);

      // TODO: 1. It seems that currentCriteria is already normalized in calculateCriteria(). we
      //          need to compare currentCriteria with 1.0 
      // TODO: 2. Current convergence criteria is not good. We need to introduce equilibrium balance
      //          criteria too along with kinematic one. 
      if (currentCriteria < convergenceCriteria) {
        converged = true;
        break;
      }
      LOG_IF(currentCriteria > 1.0e6 || std::isnan(currentCriteria), FATAL) << "The solution is diverged!";

      if (timeControl.getCurrentEquilibriumStep() >= numberOfIterations) 
        break;
    }//iterations

    LOG_IF(!converged, FATAL) << "The solution is not converged with "
        << timeControl.getCurrentEquilibriumStep() << " equilibrium iterations";
    LOG(INFO) << "Loadstep " << timeControl.getCurrentStep() << " completed with " << timeControl.getCurrentEquilibriumStep();

    // TODO: figure out why TIMED_BLOCK doesn't work here..
    // TIMED_BLOCK(t, "PostProcessor::process") {
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


double NonlinearFESolver::calculateCriteria(dVec& delta) {
  double curCriteria = 0.0;
  for (uint32 i = 0; i < delta.size(); i++) {
    curCriteria += fabs(delta[i]);
  }
  curCriteria /= delta.size();

  curCriteria = curCriteria / convergenceCriteria;
  LOG(INFO) << "Solution delta criteria = " << curCriteria;
  return curCriteria;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// LinearTransientFESolver
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
LinearTransientFESolver::LinearTransientFESolver() : FESolver() {

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

  storage->initDofs();
  setConstrainedDofs();
  storage->assignEquationNumbers();
  initSolutionData();
  

  vecR.zero();

  uint32 nEq = storage->nUnknownDofs() + storage->nMpc();
  uint32 nAll = storage->nConstrainedDofs() + storage->nUnknownDofs() + storage->nMpc();

  dVec vecUnext(nAll);
  dVec vecUnextc(vecUnext, 0, storage->nConstrainedDofs());
  dVec vecUnextsl(vecUnext, storage->nConstrainedDofs(), nEq);
  dVec vecDUnext(nAll);
  dVec vecDDUnext(nAll);

  math::SparseSymMatrix matKmod(matK->block(2)->getSparsityInfo());

  dVec rhs(nEq);

  // TODO: this hangs the program: vecU = initValue;
  for (uint32 i = 0; i < vecU.size(); i++) {
    vecU[i] = initValue;
  }

  vecDU.zero();
  vecDDU.zero();
  
  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->pre();
  }

  storage->assembleGlobalEqMatrices();
  applyBoundaryConditions(1.0);

  // dummy implementation of `matKmod = a0 * matM + a1 * matC + matK`
  for (uint32 i = 0; i < matKmod.nValues(); i++) {
    matKmod.getValuesArray()[i] = a0 * matM->block(2)->getValuesArray()[i] +
                                  a1 * matC->block(2)->getValuesArray()[i] +
                                       matK->block(2)->getValuesArray()[i];
  }

  // factorize eKmode just once as far sa we have a deal with linear system
  eqSolver->factorizeEquations(&matKmod);

  // timestepping
  while (curTime <= time1) {
    for (;;) {
      rhs = vecFsl + vecRsl;
      matBVprod(*(matM->block(2)), a0*vecUsl + a2*vecDUsl + a3*vecDDUsl, 1.0, rhs);
      matBVprod(*(matC->block(2)), a1*vecUsl + a4*vecDUsl + a5*vecDDUsl, 1.0, rhs);

      // copy constrained dofs values
      vecUnextc = vecUc;
      // solve equation system
      eqSolver->substituteEquations(&matKmod, &rhs[0], &vecUnextsl[0]);

      break;
    }//iterations

    // restore derivatives
    vecDDUnext = a0 * (vecUnext - vecU) - a2 * vecDU - a3 * vecDDU;
    vecDUnext = vecDU + a6 * vecDDUnext + a7 * vecDDU;

    vecU = vecUnext;
    vecDU = vecDUnext;
    vecDDU = vecDDUnext;

    storage->updateResults();

    LOG(INFO) << "Time " << curTime << " completed";

    for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
      postProcessors[i]->process(curTimestep);
    }

    curTimestep++;
    curTime += dt;
  } //timesteps
  LOG(INFO) << "***** SOLVED *****";

  for (size_t i = 0; i < getNumberOfPostProcessors(); i++) {
    postProcessors[i]->post(curTimestep);
  }
}

} // namespace nla3d

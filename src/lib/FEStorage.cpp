// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "FEStorage.h"
#include "elements/element.h"

namespace nla3d {
using namespace math;

FEStorage::FEStorage()  {
	numberOfNodes = 0;
	numberOfElements = 0;
	numberOfDofs = 0;

	numberOfConstrainedDofs = 0;
	numberOfUnknownDofs = 0;
	numberOfMpcEq = 0;

	material = NULL;
	KssCsT = NULL;
	Cc = NULL;
	Kcs = NULL;
	Kcc = NULL;

  elType = ElementFactory::NOT_DEFINED;

};

FEStorage::~FEStorage () {

	if (material) {
    delete material;
  }

  deleteSolutionData();
  deleteMesh();
}


void FEStorage::Kij_add(uint32 nodei, Dof::dofType dofi, uint32 nodej, Dof::dofType dofj, double value) {
	uint32 rowEq = getNodeDofEqNumber(nodei, dofi);
	uint32 colEq = getNodeDofEqNumber(nodej, dofj);
  Kij_add(rowEq, colEq, value);
}


void FEStorage::Kij_add(uint32 eqi, uint32 eqj, double value) {
  // eqi -- row equation 
  // eqj -- column equation
	if (eqi > eqj) swap(eqi, eqj);
  // As long as we have distinct blocks
  // of global matrix for constrained DoFs, MPC's lambdas and DoFs to be
  // solver we need to choose in which block (Kcc, Kcs, KssCsT)
  // the value should be added.

	if (eqi <= numberOfConstrainedDofs) {
		if (eqj <= numberOfConstrainedDofs) {
			Kcc->addValue(eqi, eqj, value);
    } else {
			Kcs->addValue(eqi, eqj - numberOfConstrainedDofs, value);
    }
  } else {
		KssCsT->addValue(eqi - numberOfConstrainedDofs, eqj - numberOfConstrainedDofs, value);
  }
}


// Cij_add is a function to add a coefficient from MPC equation to the global
// matrix.
// eq_num - number of MPC equation,
// nodej, dofj - DoFs for which the coefficient to be set.
void FEStorage::Cij_add(uint32 eq_num, uint32 nodej, Dof::dofType dofj, double coef) {
	uint32 colEq = getNodeDofEqNumber(nodej, dofj);
  Cij_add(eq_num, colEq, coef);
}


void FEStorage::Cij_add(uint32 eq_num, uint32 eqj, double coef) {
	assert(Cc);
	assert(eq_num > 0 && eq_num <= numberOfMpcEq);
	if (eqj <= numberOfConstrainedDofs) {
		Cc->addValue(eq_num, eqj, coef);
  }	else {
		//Cs
		KssCsT->addValue(eqj-numberOfConstrainedDofs, numberOfUnknownDofs+eq_num, coef);
  }
}


// Fi_add is a function to add the value to a rhs of the global system of linear
// equations. 
// NOTE: nodes index starts from 1. DoF index starts with 0.
void FEStorage::Fi_add(uint32 nodei, Dof::dofType dofi, double value) {
	uint32 rowEq = getNodeDofEqNumber(nodei, dofi);
  Fi_add(rowEq, value);
}


void FEStorage::Fi_add(uint32 eqi, double value) {
	assert(externalForces.size() > 0);
	assert(eqi <= numberOfDofs);
	externalForces[eqi - 1] += value;
}


// zeroK function is used to set values in the global stiffnes matrix to zero.
void FEStorage::zeroK() {
	assert(KssCsT);
	KssCsT->zero();
	Kcc->zero();
	Kcs->zero();
  if (Cc) {
    Cc->zero();
  }
}

// zeroF:
// set zeros to a rhs vector of a global eq. system.
void FEStorage::zeroF() {
	assert(externalForces.size() > 0);

  std::fill(externalForces.begin(),externalForces.end(), 0.0); 
  std::fill(solutionDeltaValues.begin(), solutionDeltaValues.end(), 0.0); 

  if (numberOfMpcEq > 0) {
    std::fill(mpcConstantValues.begin(),mpcConstantValues.end(), 0.0); 
  }
}

// getGlobalEqMatrix:
math::SparseSymmetricMatrix* FEStorage::getGlobalEqMatrix() {
	assert(KssCsT && Kcc && Kcs);
	return KssCsT;
}

// prepare and pass back a rhs of the global eq. system.
// As we exclude constrained DoFs from a eq. system we need to modify rhs of the
// system.
double* FEStorage::getGlobalEqRhs () {
	assert(rhs.size() > 0);
	double *KcsTdqc = new double[numberOfUnknownDofs];
	Kcs->transpose_mult_vec(&constrainedDofDeltaValues[0], KcsTdqc);
	for (uint32 i=0; i < numberOfUnknownDofs; i++) {
    // was:
    //Kcs->transpose_mult_vec_i(&constrainedDofDeltaValues[0],i+1);
		unknownDofRhs[i] = unknownDofExternalForces[i] - KcsTdqc[i];
  }
	for (uint32 i=0; i < numberOfMpcEq; i++) {
		mpcEquationRhs[i] = mpcConstantValues[i] - Cc->mult_vec_i(&constrainedDofDeltaValues[0], i + 1);
  }
	delete[] KcsTdqc;
	return &rhs[0];
}

double* FEStorage::getGlobalEqUnknowns () {
	assert(solutionDeltaValues.size() > 0);
	return &unknownDofDeltaValues[0];
}

void FEStorage::assembleGlobalEqMatrix() {
  TIMED_SCOPE(t, "assembleGlobalEqMatrix");
  LOG(INFO) << "Start formulation of elements stiffness matrices ( " << getNumberOfElements() << " elements)";
	CHECK(KssCsT && Kcc && Kcs);

  if (!KssCsT->isInTrainingMode()) {
    zeroK();
    zeroF();
  }

  for (uint32 el = 0; el < getNumberOfElements(); el++) {
    elements[el]->build();
  }
  //t.checkpoint("Element::build()");

  // because of non-linear MPC we need to update
  // MPC coefficients every step
  //
  for (size_t i = 0; i < mpcCollections.size(); i++) {
    mpcCollections[i]->update();
    //mpcCollections[i]->printEquations(std::cout);
  }
  //t.checkpoint("MpcCollection::update()");

  uint32 eq_num = 1;
  list<Mpc*>::iterator mpc = mpcs.begin();
  while (mpc != mpcs.end()) {
    mpcConstantValues[eq_num-1] = (*mpc)->b;
    list<MpcTerm>::iterator token = (*mpc)->eq.begin();
    while (token != (*mpc)->eq.end()) {
      // TODO: now we support only MPC to nodal DoFs..
      Cij_add(eq_num, token->node, token->node_dof, token->coef);
      token++;
    }
    eq_num++;
    mpc++;
  }

  // finishSparseMatricesTraining
  // for details on it see a comment with startSparseMatricesTraining mark
  if (KssCsT->isInTrainingMode()) {
    TIMED_SCOPE(t, "SparseMatrix::stopTraining");
    assert(Kcc->isInTrainingMode() && Kcs->isInTrainingMode());
    if (numberOfMpcEq) {
      assert(Cc->isInTrainingMode());
      Cc->stopTraining();
    }
    KssCsT->stopTraining();
    
    Kcc->stopTraining();
    Kcs->stopTraining();
  }
}

uint32 FEStorage::getNumberOfNodes () {
  return static_cast<uint32> (nodes.size());
}

uint32 FEStorage::getNumberOfElements () {
  return static_cast<uint32> (elements.size());
}

uint32 FEStorage::getNumberOfUnknownDofs() {
  return numberOfUnknownDofs;
}

uint32 FEStorage::getNumberOfConstrainedDofs() {
  return static_cast<uint32> (constraints.size());
}

uint32 FEStorage::getNumberOfMpcEq() {
  return static_cast<uint32> (mpcs.size());
}

void FEStorage::addNodeDof(uint32 node, std::initializer_list<Dof::dofType> _dofs) {
  assert(nodeDofs.getNumberOfEntities() > 0);
  nodeDofs.addDof(node, _dofs);
}

void FEStorage::addElementDof(uint32 el, std::initializer_list<Dof::dofType> _dofs) {
  assert(elementDofs.getNumberOfEntities() > 0);
  elementDofs.addDof(el, _dofs);
}


bool FEStorage::isElementDofUsed (uint32 el, Dof::dofType dof) {
  assert(elementDofs.getNumberOfEntities() > 0);
  return elementDofs.isDofUsed(el, dof);
}


bool FEStorage::isNodeDofUsed (uint32 node, Dof::dofType dof) {
  assert(nodeDofs.getNumberOfEntities() > 0);
  return nodeDofs.isDofUsed(node, dof);
}


uint32 FEStorage::getElementDofEqNumber(uint32 el, Dof::dofType dof) {
  return getElementDof(el, dof)->eqNumber;
}


uint32 FEStorage::getNodeDofEqNumber(uint32 node, Dof::dofType dof) {
  return getNodeDof(node, dof)->eqNumber;
}


Dof* FEStorage::getElementDof (uint32 el, Dof::dofType dof) {
  assert(elementDofs.getNumberOfUsedDofs() > 0);
  return elementDofs.getDof(el, dof);
}

Dof* FEStorage::getNodeDof (uint32 node, Dof::dofType dof) {
  assert(nodeDofs.getNumberOfUsedDofs() > 0);
  return nodeDofs.getDof(node, dof);
}


double FEStorage::getNodeDofSolution (uint32 node, Dof::dofType dof) {
  assert(solutionValues.size() > 0);
  return solutionValues[getNodeDofEqNumber(node, dof) - 1];
}


double FEStorage::getElementDofSolution (uint32 el, Dof::dofType dof) {
  assert(solutionValues.size() > 0);
  return solutionValues[getElementDofEqNumber(el, dof) - 1];
}


double FEStorage::getNodeDofSolutionDelta (uint32 node, Dof::dofType dof) {
  assert(solutionDeltaValues.size() > 0);
  return solutionDeltaValues[getNodeDofEqNumber(node, dof) - 1];
}


double FEStorage::getElementDofSolutionDelta (uint32 el, Dof::dofType dof) {
  assert(solutionDeltaValues.size() > 0);
  return solutionDeltaValues[getElementDofEqNumber(el, dof) - 1];
}


double FEStorage::getReaction(uint32 node, Dof::dofType dof) {
  if (isNodeDofUsed(node, dof)) {
    uint32 eq_num = getNodeDofEqNumber(node,dof);
    return getReaction(eq_num);
  } else {
    return 0.0;
  }
}


double FEStorage::getReaction(uint32 eq) {
  assert(constrainedDofReactions.size() > 0);
  //assert(eq_num > 0 && eq_num <= numberOfConstrainedDofs);
  if (eq > 0 && eq <= numberOfConstrainedDofs) {
    return constrainedDofReactions[eq-1];
  } else {
    return 0.0;
  }
}

Material* FEStorage::getMaterial () {
	assert(material);
	return material;
}

// _nn > 0 (numbers start from 1)
Node& FEStorage::getNode(uint32 _nn) {
	assert(_nn <= numberOfNodes);
	return *nodes[_nn-1];
}

void FEStorage::getElementNodes(uint32 el, Node** node_ptr)
{
	assert(el <= numberOfElements);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		node_ptr[i] = nodes[elements[el-1]->getNodeNumber(i)-1];
}

// n starts from 1
// prt array always has 3 elements
void FEStorage::getNodePosition(uint32 n, double* ptr, bool deformed)
{
	assert(n > 0 && n <= numberOfNodes);
	for (uint16 i = 0; i < 3; i++) {
		ptr[i] = nodes[n-1]->pos[i];
  }
	if (deformed) {
    if (isNodeDofUsed(n, Dof::UX)) {
      ptr[0] += getNodeDofSolution(n, Dof::UX);
    }

    if (isNodeDofUsed(n, Dof::UY)) {
      ptr[1] += getNodeDofSolution(n, Dof::UY);
    }

    if (isNodeDofUsed(n, Dof::UZ)) {
      ptr[2] += getNodeDofSolution(n, Dof::UZ);
    }
	}
}

// _en > 0 (numbers start from 1)
Element& FEStorage::getElement(uint32 _en) {
	assert(_en <= numberOfElements);
	return *(elements[_en-1]);
}

FEComponent* FEStorage::getFEComponent(size_t i) {
  assert(i < feComponents.size());
  return feComponents[i];
}

FEComponent* FEStorage::getFEComponent(const std::string& name) {
  for (size_t i = 0; i < feComponents.size(); i++) {
    if (name.compare(feComponents[i]->name) == 0) {
      return feComponents[i];
    }
  }
  LOG(WARNING) << "Can't find a component with name " <<  name;
  return NULL;
}

void FEStorage::addFixation (int32 n, Dof::dofType dof, const double value) {
  // FEStorage believes that all BC_dof_constraint are distinct!
  // Do not pass to it the same BC twice!
  constraints.push_back(BC_dof_constraint(n, dof, value));
}

void FEStorage::addForce(int32 n, Dof::dofType dof, const double value) {
  forces.push_back(BC_dof_force (n, dof, value));
}

void FEStorage::addBoundaryCondition (BC_dof_constraint &bc) {
  // FEStorage believes that all BC_dof_constraint are distinct!
  // Do not pass to it the same BC twice!
	constraints.push_back(bc);
}

void FEStorage::addBoundaryCondition (BC_dof_force &bc) {
  // FEStorage believes that all BC_dof_force are distinct!
  // Do not pass to it the same BC twice!
	forces.push_back(bc);
}

void FEStorage::addMpc (Mpc* mpc) {
  assert (mpc);
  assert (mpc->eq.size() > 0);
	mpcs.push_back(mpc);
}

void FEStorage::addMpcCollection (MpcCollection* mpcCol) {
  assert(mpcCol);
  mpcCollections.push_back(mpcCol);
}

void FEStorage::addFEComponent (FEComponent* comp) {
	assert(comp);
	feComponents.push_back(comp);
}

void FEStorage::addNode (Node* node) {
	numberOfNodes++;
  //Node() fires Vec<3> constructor, thus Node coordinates are (0,0,0) by default
  //TODO: try-catch of memory overflow
	nodes.push_back(node);
}

void FEStorage::createNodes (uint32 _nn) {
	deleteNodes();
	numberOfNodes = _nn;
  //Node() fires Vec<3> constructor, thus Node coordinates are (0,0,0) by default
  //TODO: try-catch of memory overflow
  nodes.reserve(_nn);
  for (uint32 i = 0; i < _nn; i++) {
    nodes.push_back(new Node);
  }
}

void FEStorage::addElement (Element* el) {
	numberOfElements++;
  //TODO: try-catch of memory overflow
	elements.push_back(el);
  Element::storage = this;
}

//createElements(_en)
void FEStorage::createElements(uint32 _en) {
  deleteElements();
	numberOfElements = _en;
  //TODO: catch if not enough memory
  elements.reserve(_en);
  ElementFactory::createElements (elType, numberOfElements, elements); 
  Element::storage = this;
  for (uint32 i = 0; i < _en; i++) {
    //access elNum protected values as friend
    elements[i]->elNum = i+1;
  }
}

void FEStorage::deleteMesh () {
	deleteElements();
	deleteNodes();
	constraints.clear();
	forces.clear();
  deleteMpcs();
  deleteMpcCollections();
  deleteFeComponents();
	numberOfNodes = 0;
}

void FEStorage::deleteNodes () {
  for (uint32 i = 0; i < numberOfNodes; i++) {
    delete nodes[i];
  }
  nodes.clear();
  numberOfNodes = 0;
}
// delete element table
void FEStorage::deleteElements() {
  for (uint32 i = 0; i < numberOfElements; i++) {
    delete elements[i];
  }
  elements.clear();
  numberOfElements = 0;
}

void FEStorage::deleteMpcs() {
  list<Mpc*>::iterator mpc = mpcs.begin();
  while (mpc != mpcs.end()) {
    delete *mpc;
    mpc++;
  }
  numberOfMpcEq = 0;
}

void FEStorage::deleteMpcCollections () {

  for (size_t i = 0; i < mpcCollections.size(); i++) {
    delete mpcCollections[i];
  }
  mpcCollections.clear();
}

void FEStorage::deleteFeComponents () {
  for (size_t i = 0; i < feComponents.size(); i++) {
    delete feComponents[i];
  }
  feComponents.clear();
}

void FEStorage::deleteDofArrays() {
  elementDofs.clearDofTable();
  nodeDofs.clearDofTable();
}

void FEStorage::deleteSolutionData() {
	numberOfDofs = 0;

  numberOfConstrainedDofs = 0;
	numberOfMpcEq = 0;

	numberOfUnknownDofs = 0;

  if (Cc) {
    delete Cc;
    Cc = NULL;
  }
  if (KssCsT) {
    delete KssCsT;
    KssCsT = NULL;
  }
  if (Kcs) {
    delete Kcs;
    Kcs = NULL;
  }
  if (Kcc) {
    delete Kcc;
    Kcc = NULL;
  }
  
  deleteDofArrays();

  solutionValues.clear();
  dofValues = solutionValues.end();
  constrainedDofValues = solutionValues.end();
  unknownDofValues = solutionValues.end();
  mpcLagrangianValues = solutionValues.end();
	
  solutionDeltaValues.clear();
	dofDeltaValues = solutionDeltaValues.end();
	constrainedDofDeltaValues = solutionDeltaValues.end();
  unknownDofDeltaValues = solutionDeltaValues.end();
  mpcLagrangianDeltaValues = solutionDeltaValues.end();
	
	constrainedDofReactions.clear();

	externalForces.clear();
	constrainedDofExternalForces = externalForces.end();
	unknownDofExternalForces = externalForces.end();

  mpcConstantValues.clear();

	rhs.clear();
	unknownDofRhs = rhs.end();
	mpcEquationRhs = rhs.end();

	numberOfDofs = 0;
	numberOfConstrainedDofs = 0;
	numberOfUnknownDofs = 0;
}

void FEStorage::listFEComponents () {
  for (size_t i = 0; i < feComponents.size(); i++) {
    LOG(INFO) << *feComponents[i];
  }
}

bool FEStorage::initializeSolutionData () {
  TIMED_SCOPE(t, "initializeSolutionData");
	if (getNumberOfNodes() < 1 && getNumberOfElements() < 1) {
		LOG(ERROR) << "numberOfNodes or numberOfElements is zero";
    exit(1);
	}
  
  nodeDofs.initDofTable(getNumberOfNodes());
  elementDofs.initDofTable(getNumberOfElements());

	for (uint32 el = 0; el < getNumberOfElements(); el++) {
    elements[el]->pre();
  }
  //t.checkpoint("Element::pre()");

  for (size_t i = 0; i < mpcCollections.size(); i++) {
    mpcCollections[i]->pre();
    mpcCollections[i]->registerMpcsInStorage();
  }
  //t.checkpoint("MpcCollection::pre()");

  // Total number of dofs (only registered by elements)
	numberOfDofs = elementDofs.getNumberOfUsedDofs() + nodeDofs.getNumberOfUsedDofs();

  //TODO: make it possible to work without any constrained dofs
  // (in case of only MPC)
  numberOfConstrainedDofs = static_cast<uint32> (constraints.size());
	numberOfMpcEq = static_cast<uint32> (mpcs.size());

  // numberOfUnknownDofs - number of Dof need to be found on every step
	numberOfUnknownDofs = numberOfDofs - numberOfConstrainedDofs;
  // In nla3d solution procedure there are 3 distinguish types of unknowns. First one "c" - constrained degress of freedom,
  // and consequently known at solution time. Second one "s" - degrees of freedom need to be found (solved).
  // And last one "lambda" - lagrange multipliers for applied MPC constraints.
  // Cc [numberOfMpcEq x numberOfConstrainedDofs]  - part of a global stiffness matrix with MPC coefficients for constrained dofs
  // A global System of Linera Algebraic Equations:
  //
  //  |  Kcc  |  Kcs  | Cc^T |   |   qc  |   |constrainedDofExternalForces|
  //  |-------|--------------|   |-------|   |      |
  //  | Kcs^T |              | * |   qs  | = |unknownDofExternalForces|
  //  |-------|    KssCsT    |   |-------|   |      |
  //  |  Cc   |              |   | lambda|   |mpcConstantValues |
  //  
  //  1. But really only
  //
  //  KssCsT * [qs; lambda]^T = [unknownDofRhs; mpcEquationRhs]^T
  //
  //  to be solved by eq. solver. 
  //
  //  where:
  //  unknownDofRhs = unknownDofExternalForces - Kcs^T*qc
  //  mpcEquationRhs = mpcConstantValues - Cc*qc
  //
  //  2. Then to restore reaction forces:
  //  constrainedDofExternalForces = Kcc*qc + Kcs*qs + Cc^T*lambda
  //
	if (numberOfMpcEq) {
    Cc = new math::SparseMatrix(numberOfMpcEq, numberOfConstrainedDofs);
  }
	KssCsT = new math::SparseSymmetricMatrix(numberOfUnknownDofs + numberOfMpcEq);
	Kcs = new math::SparseMatrix(numberOfConstrainedDofs, numberOfUnknownDofs);
	Kcc = new math::SparseSymmetricMatrix(numberOfConstrainedDofs);

	// Fill elementsDofs and nodeDofs with isConstrained information
  // and then give them an equation numbers
	uint32 next_eq_solve = numberOfConstrainedDofs+1;
	uint32 next_eq_const = 1;
	list<BC_dof_constraint>::iterator p = constraints.begin();
	while (p != constraints.end()) {
    //TODO: currently we work only with nodal DoF constraints, but need to fix it!
    getNodeDof(p->node, p->node_dof)->isConstrained = true;
		p++;
	}

  for (uint32 i = 1; i <= getNumberOfElements(); i++) {
    for (uint16 it = 0; it < Dof::numberOfDofTypes; it++) {
      Dof::dofType t = static_cast<Dof::dofType> (it);
      Dof* d = elementDofs.getDof(i, t);
      if (d) {
        if (d->isConstrained) {
          d->eqNumber = next_eq_const++;
        } else {
          d->eqNumber = next_eq_solve++;
        }
      }
    }
  }


  for (uint32 i = 1; i <= getNumberOfNodes(); i++) {
    for (uint16 it = 0; it < Dof::numberOfDofTypes; it++) {
      Dof::dofType t = static_cast<Dof::dofType> (it);
      Dof* d = nodeDofs.getDof(i, t);
      if (d) {
        if (d->isConstrained) {
          d->eqNumber = next_eq_const++;
        } else {
          d->eqNumber = next_eq_solve++;
        }
      }
    }
  }

  assert(next_eq_const - 1 == numberOfConstrainedDofs);
  assert(next_eq_solve - 1 == numberOfDofs);

	// initialize solutionValues vector and setup pointers to specific parts of the vector
  solutionValues.assign(numberOfDofs+numberOfMpcEq, 0.0);
  dofValues = solutionValues.begin();
  constrainedDofValues = solutionValues.begin();
  unknownDofValues = solutionValues.begin() + numberOfConstrainedDofs;
  mpcLagrangianValues = solutionValues.begin() + numberOfDofs;
	
  solutionDeltaValues.assign(numberOfDofs+numberOfMpcEq, 0.0);

	dofDeltaValues = solutionDeltaValues.begin();
	constrainedDofDeltaValues = solutionDeltaValues.begin();
  unknownDofDeltaValues = solutionDeltaValues.begin() + numberOfConstrainedDofs;
  mpcLagrangianDeltaValues = solutionDeltaValues.begin() + numberOfDofs;
	
	constrainedDofReactions.assign(numberOfConstrainedDofs, 0.0);

	externalForces.assign(numberOfDofs, 0.0);
	constrainedDofExternalForces = externalForces.begin();
	unknownDofExternalForces = externalForces.begin() + numberOfConstrainedDofs;


	if (numberOfMpcEq) {
		mpcConstantValues.assign(numberOfMpcEq, 0.0);
	}

	// size of rhs [numberOfUnknownDofs + numberOfMpcEq]
	rhs.assign(numberOfUnknownDofs + numberOfMpcEq, 0.0);
	unknownDofRhs = rhs.begin();
	mpcEquationRhs = rhs.begin() + numberOfUnknownDofs;

  std::stringstream ss;
  ss << " Number of nodal DoFs: " << nodeDofs.getNumberOfUsedDofs();
  ss << std::endl;
  ss << " Number of element DoFs: " << elementDofs.getNumberOfUsedDofs();

  LOG(INFO) << ss.str();
	LOG(INFO) << "DoFs = " << numberOfDofs << ", constrained DoFs = " <<  numberOfConstrainedDofs << ", MPC eq. = "
      <<  numberOfMpcEq << ", TOTAL eq. = " << numberOfUnknownDofs + numberOfMpcEq;

  //t.checkpoint("Prepearing FEStorage data for solution");
	if (!material) {
		LOG(WARNING) << "FEStorage::initializeSolutionData: material isn't defined";
	}

  // startSparseMatricesTraining
  // Because apriori FEStorage doesn't know about which elements in sparse matrices are non zero.
  // To estimate this here is a specual procedure named Training (see SparseMatrix classes).
  // Training need to be done on the first global matrix assembling. After training is done a data
  // obout sparse matrix elements are astored in special arrays in SparseMatrix classes. After this
  // an action of writing to zero (non-existent) element into a SparseMatrix will be enede with errors.
  // That means that elements and MPC can't change number of used DoFs during solving procedures. 
  // That is a bad thing in, for example, non-linear MPC, where number of DoFs involved could be differ
  // during the solution.
  // To see where is training procedures ended look for finishSparseMatricesTraining comment.
	if (numberOfMpcEq) {
    Cc->startTraining();
  }
	KssCsT->startTraining();

	Kcc->startTraining();
	Kcs->startTraining();
	return true;
}

void FEStorage::updateSolutionResults() {
  TIMED_SCOPE(t, "updateSolutionResults");
  // add delta to DoFs values from current solved iteration
	for (size_t i = 0; i < numberOfDofs; i++) {
		dofValues[i] += dofDeltaValues[i];
  }

  // lagrangian from MPC are updated from current solved step
	for (uint32 i = 0; i < numberOfMpcEq; i++) {
    mpcLagrangianValues[i] = mpcLagrangianDeltaValues[i];
  }
	
  // also update reactions for constrained DoFs
	for (uint32 i = 0; i < numberOfConstrainedDofs; i++) {
		constrainedDofReactions[i] = Kcs->mult_vec_i(&unknownDofDeltaValues[0], i + 1) +
        Kcc->mult_vec_i(&constrainedDofDeltaValues[0], i + 1) - constrainedDofExternalForces[i];

    // in case of constrained DoFs are involved in MPC equations it's needed to add reactions from the MPCs
		if (numberOfMpcEq && Cc->getNumberOfValues()) {
			constrainedDofReactions[i] += Cc->transpose_mult_vec_i(&mpcLagrangianDeltaValues[0], i + 1);
    }
	}
  //t.checkpoint("Sum solutions. Finding reactions");

  // calculate element's update procedures (calculate stresses, strains, ..)
  for (uint32 el = 0; el < getNumberOfElements(); el++) {
    elements[el]->update();
  }
  //t.checkpoint("Element::update");
}


void FEStorage::applyBoundaryConditions (double time, double timeDelta) {
  TIMED_SCOPE(t, "applyBoundaryConditions");
  LOG(INFO) << "Aplaying boundary conditions.. (" << forces.size() << " nodal loads and "
       << constraints.size() << "nodal fixations)";	
	// fill nodal forces
	list<BC_dof_force>::iterator bc_force = forces.begin();
	while (bc_force != forces.end())	{
		Fi_add(bc_force->node, bc_force->node_dof, bc_force->value*time);
		bc_force++;
	}
  //t.checkpoint("apply forces");

	// fill nodal displacements (kinematic constraints)
	list<BC_dof_constraint>::iterator bc_dof = constraints.begin();
	while (bc_dof != constraints.end()) {
    //TODO: now support only nodal dofs..
		uint32 eq_num = getNodeDofEqNumber(bc_dof->node, bc_dof->node_dof);
    // To be sure that constrained DoF lays in numberOfConstrainedDofs part
    assert (eq_num - 1 < numberOfConstrainedDofs);
		dofDeltaValues[eq_num-1] = bc_dof->value*timeDelta;
		bc_dof++;
	}
  //t.checkpoint("apply fixations");
}

bool readCdbFile(const char *filename, FEStorage *storage)
{
	uint32 n_number, en;
	ifstream file(filename);
	if (!file) {
		LOG(WARNING) << "Can't open input file " << filename;
		return false;
	}
	storage->deleteMesh();
	char buf[1024]="";
	while (file.getline(buf, 1024))
	{
		vector<string> vec = read_tokens(buf);
		if (vec[0].compare("NBLOCK") == 0)
		{
    //NBLOCK,6,SOLID,     9355,     9355
    //(3i9,6e20.13)
    //        1        0        0 7.0785325971794E+01 6.5691449317818E+01-3.6714639015390E+01
			uint32 max_n_number = atoi(vec[6].c_str());
			n_number= atoi(vec[8].c_str());
      if (max_n_number != n_number) {
        LOG(ERROR) << "NBLOCK: maximum node number is " << max_n_number
            << "but number of nodes is " << n_number << ". Note that nla3d needs compressed numbering "
            << "for nodes and elements";
        exit(1);
      }
			storage->createNodes(n_number);
			file.getline(buf, 1024);
      string buf_str(buf);
      // we need to take a format of columns "3i9"
      size_t start = buf_str.find("i")+1;
      size_t  stop = buf_str.find(",");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str());

			for (uint32 i=1; i<=n_number; i++) {
				file.getline(buf, 1024);
				size_t len=strlen(buf);
				for (uint16 j=0; j<3;j++)
					if (len>=3*frmt+20*(j+1))
            //note that last column in NBLOCK table could be avoided if Z=0.0
            //but storage->createNodes(n_number) initialize the node table with (0,0,0)
						storage->getNode(i).pos[j] = atof(string((char*) (buf+3*frmt+20*j),20).c_str());
			}
		}//NBLOCK
		else if (vec[0].find("EBLOCK")!=vec[0].npos)
		{  
      //EBLOCK,19,SOLID,      7024,      7024
      //(19i9)
			en = atoi(vec[6].c_str());
      if (en != atoi(vec[8].c_str())) {
        LOG(ERROR) << "EBLOCK: maximum element number is " << en << ", but number of element "
            << "is different. Note that nla3d needs compressed numbering for nodes and elements";
        exit(1);
      }
			storage->createElements(en);
			file.getline(buf, 1024);
      // we need to take a format of columns "3i9"
      // in Ansys 12 here is 8 symbols per number (19i8), but in ansys 15 (19i9) is used. 
      string buf_str(buf);
      size_t start = buf_str.find("i")+1;
      size_t stop = buf_str.find(")");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str()); 
			for (uint32 i=1; i<=en; i++)
			{
				file.getline(buf, 1024);
				size_t len=strlen(buf);
        // TODO: here is a problem.. We need to work good with both dos and unix endings
//        if (buf[len-1] == '\r') {
//          len--;
//          buf[len-1] = 0;
//        }
        //TODO: It seems that getline keeps windows line ending
        if (len != 11*frmt+frmt*Element::n_nodes()) {
          LOG_N_TIMES(10, WARNING) << "in EBLOCK for element " << i << " the number of nodes provided is not equal to " << Element::n_nodes();
        }
				for (uint16 j=0; j<Element::n_nodes();j++)
					if (len>=11*frmt+frmt*(j+1))
						storage->getElement(i).getNodeNumber(j) = atoi(string((char*) (buf+11*frmt+frmt*j),frmt).c_str());
			}
		}//EBLOCK
		else if (vec[0].find('D')!=vec[0].npos && vec[0].length()==1)
		{
				BC_dof_constraint bnd;
				bnd.node = atoi(vec[2].c_str());
				bnd.value = atof(vec[6].c_str());
        bnd.node_dof = Dof::label2dofType(vec[4]);
        storage->addBoundaryCondition(bnd);
		}//D
    else if (vec[0].compare("CE") == 0)
    {
      //How MPC looks like this in cdb file:
      //CE,R5.0,DEFI,       2,       1,  0.00000000    
      //CE,R5.0,NODE,      1700,UX  ,  1.00000000    ,      1700,UZ  ,  1.00000000  
      Mpc* mpc = new Mpc;
      mpc->b = atoi(vec[10].c_str()); //rhs of MPC equation
      uint16 n_terms = atoi(vec[6].c_str()); //number of terms in equation
      //debug("MPC link: %d terms, b = %f", n_terms, mpc.b);
      while (n_terms > 0)
      {
        file.getline(buf, 1024);
        vector<string> vec = read_tokens(buf);
        uint16 place = 6;
        for (int i=0; i < max((uint16) n_terms, (uint16) 2); i++) 
        {
          uint32 node = atoi(vec[place].c_str());
          Dof::dofType dof = Dof::label2dofType(vec[place+2]);
          double coef = atof(vec[place+4].c_str());
          //debug("%d term: node %d, dof %d, coef = %f", i, node, dof, coef);
          mpc->eq.push_back(MpcTerm(node,dof,coef));
          place += 6;
          n_terms--;
        }
      }
			storage->addMpc(mpc);
    }//CE (MPC)
    else if (vec[0].compare("CMBLOCK") == 0) {
       if (vec[2].c_str()[0] != '_') {
         FEComponent* comp = new FEComponent();
         comp->name = vec[2];
      //CMBLOCK,BOTTOM_SIDE,NODE,      17  ! users node component definition
      //(8i10)
      //      5037     -5341      6330     -6352      6355      6357      6433     -6456
      //      6459      6470     -6473      6537     -6556      6566     -6569      6633
      //     -6652
        comp->type = FEComponent::typeFromString(vec[4]);

        size_t numRanges = atoi(vec[6].c_str());
        file.getline(buf, 1024);
        // we need to take a format of columns "8i10"
        // read how many records in a single row
        string buf_str(buf);
        size_t start = buf_str.find("(")+1;
        size_t stop = buf_str.find("i");
        string tmp;
        tmp.assign((char*) (buf + start), stop-start);
        uint16 numRangesInRow = atoi(tmp.c_str()); 
        // read how many symbols dedicated to a number
        start = buf_str.find("i")+1;
        stop = buf_str.find(")");
        tmp.assign((char*) (buf + start), stop-start);
        uint16 frmt = atoi(tmp.c_str()); 
        uint16 in_row = 0;
        uint16 all = 0;

        file.getline(buf, 1024);
        vector<int32> rangesVec;
        rangesVec.reserve(numRanges);
        size_t numEntity = 0;
        while (all < numRanges) {
           if (in_row == numRangesInRow) {
              file.getline(buf, 1024);
              in_row = 0;
           }
           tmp.assign((char*) (buf+in_row*frmt),frmt);
           rangesVec.push_back(atoi(tmp.c_str())); 
           if (rangesVec[all] > 0) {
              numEntity++;
           } else {
             //4 -9 : 4 5 6 7 8 9
              assert(all > 0);
              assert(rangesVec[all-1] > 0);
              numEntity += -rangesVec[all] - rangesVec[all-1];
           }
           in_row ++;
           all ++;
        }
        comp->list.reserve(numEntity);
        all = 0;
        while (all < numRanges) {
          if (rangesVec[all] > 0) {
            comp->list.push_back(rangesVec[all]);
          } else {
            for (uint32 i = rangesVec[all-1] + 1; i < static_cast<uint32> (-rangesVec[all]+1); i++) {
              //TODO: need assert for overflow
              comp->list.push_back(i);
            }
          }
          all++;
        }
        storage->addFEComponent(comp);
      }
    } //CMBLOCK
    //TODO: add FX FY FZ
	}
	file.close();
	return true;
}





} // namespace nla3d 

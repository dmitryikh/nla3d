// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "FEStorage.h"
#include "elements/element.h"

namespace nla3d {
using namespace math;


FEStorage::FEStorage()  {
};


FEStorage::~FEStorage () {
	if (material) {
    delete material;
  }

  deleteSolutionData();
  deleteMesh();
}


void FEStorage::assembleGlobalEqMatrices() {
	assert(matK);
  assert(matK->isCompressed());

  TIMED_SCOPE(t, "assembleGlobalEqMatrix");
  LOG(INFO) << "Start formulation of global eq. matrices ( " << nElements() << " elements)";

  zeroK();
  zeroF();

  for (uint32 el = 0; el < nElements(); el++) {
    elements[el]->buildK();
  }
  //t.checkpoint("Element::build()");

  // because of non-linear MPC we need to update
  // MPC coefficients every step
  for (size_t i = 0; i < mpcCollections.size(); i++) {
    mpcCollections[i]->update();
    // DEBUG:
    // mpcCollections[i]->printEquations(std::cout);
  }
  //t.checkpoint("MpcCollection::update()");

  // loop over mpc equations and add corresponding terms into global eq. system
  for (auto& mpc : mpcs) {
    uint32 eq_num = mpc->eqNum;
    assert(eq_num > 0);
    assert(eq_num <= vecF.size());
    assert(eq_num <= nDofs() + nMpc());
    vecF[eq_num - 1] = mpc->b;

    for (auto& term : mpc->eq) {
      addValueMPC(mpc->eqNum, term.node, term.node_dof, term.coef);
    }
  }

  // if it was demanded to have transient matrices (C and M)
  if (transient) {
    assert(matC && matM);
    assert(matC->isCompressed());
    assert(matM->isCompressed());

    zeroC();
    zeroM();

    for (uint32 el = 0; el < nElements(); el++) {
      elements[el]->buildC();
      elements[el]->buildM();
    }
  }

}


void FEStorage::setConstrainedNodeDof(uint32 node, Dof::dofType dtype) {
  Dof* dof = nodeDofs.getDof(node, dtype);
  assert(dof);
  if (dof->isConstrained) {
    return;
  } else {
    dof->isConstrained = true;
    _nConstrainedDofs++;
  }
}


void FEStorage::setConstrainedElementDof(uint32 el, Dof::dofType dtype) {
  Dof* dof = elementDofs.getDof(el, dtype);
  assert(dof);
  if (dof->isConstrained) {
    return;
  } else {
    dof->isConstrained = true;
    _nConstrainedDofs++;
  }
}


double FEStorage::getReaction(uint32 node, Dof::dofType dof) {
  if (isNodeDofUsed(node, dof)) {
    uint32 eq_num = getNodeDofEqNumber(node, dof);
    return getReaction(eq_num);
  } else {
    return 0.0;
  }
}


double FEStorage::getReaction(uint32 eq) {
  assert(eq > 0 && eq <= vecR.size());
  // NOTE: vecR contains not only reactions but also external forces (for eq > nConstrainedDofs())
  // so this function will return external force for eq > nConstrainedDofs()
  return vecR[eq - 1];
}

Material* FEStorage::getMaterial() {
	assert(material);
	return material;
}


void FEStorage::getElementNodes(uint32 el, Node** node_ptr) {
	assert(el <= nElements());
  Element* elp = elements[el-1];
	for (uint16 i=0; i<elp->getNNodes(); i++)
		node_ptr[i] = nodes[elp->getNodeNumber(i)-1];
}


// prt array always has 3 elements
void FEStorage::getNodePosition(uint32 n, double* ptr, bool deformed) {
	assert(n > 0 && n <= nNodes());
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
  //Node() fires Vec<3> constructor, thus Node coordinates are (0,0,0) by default
  //TODO: try-catch of memory overflow
	nodes.push_back(node);
}


std::vector<uint32> FEStorage::createNodes (uint32 _nn) {
  //Node() fires Vec<3> constructor, thus Node coordinates are (0,0,0) by default
  //TODO: try-catch of memory overflow
  std::vector<uint32> newIndexes;
  newIndexes.reserve(_nn);
  uint32 nextNumber = nodes.size() + 1;

  nodes.reserve(elements.size() + _nn);

  for (uint32 i = 0; i < _nn; i++) {
    Node* pnode = new Node;
    nodes.push_back(pnode);
  }

  for (uint32 i = nextNumber; i <= nodes.size(); i++) {
    newIndexes.push_back(i);
  }
  return newIndexes;
}


void FEStorage::addElement (Element* el) {
  el->storage = this;
  el->elNum = nElements() + 1;
  //TODO: try-catch of memory overflow
	elements.push_back(el);
}


std::vector<uint32> FEStorage::createElements(uint32 _en, ElementType elType) {
  //TODO: catch if not enough memory
  std::vector<uint32> newIndexes;
  newIndexes.reserve(_en);
  uint32 nextNumber = elements.size() + 1;
  ElementFactory::createElements(elType, _en, elements); 
  for (uint32 i = nextNumber; i <= elements.size(); i++) {
    //access elNum protected values as friend
    elements[i - 1]->elNum = i;
    elements[i - 1]->storage = this;
    newIndexes.push_back(i);
  }
  return newIndexes;
}


void FEStorage::deleteMesh() {
	deleteElements();
	deleteNodes();
  deleteMpcs();
  deleteMpcCollections();
  deleteFeComponents();
  topology.clear();
}


void FEStorage::deleteNodes() {
  for (uint32 i = 0; i < nNodes(); i++) {
    delete nodes[i];
  }
  nodes.clear();
}


void FEStorage::deleteElements() {
  for (uint32 i = 0; i < nElements(); i++) {
    delete elements[i];
  }
  elements.clear();
}


void FEStorage::deleteMpcs() {
  list<Mpc*>::iterator mpc = mpcs.begin();
  while (mpc != mpcs.end()) {
    delete *mpc;
    mpc++;
  }
}


void FEStorage::deleteMpcCollections() {

  for (size_t i = 0; i < mpcCollections.size(); i++) {
    delete mpcCollections[i];
  }
  mpcCollections.clear();
}


void FEStorage::deleteFeComponents() {
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
	_nDofs = 0;
  _nConstrainedDofs = 0;
	_nUnknownDofs = 0;

  if (matK) {
    delete matK;
    matK = nullptr;
  }
  if (matC) {
    delete matC;
    matC = nullptr;
  }
  if (matM) {
    delete matM;
    matM = nullptr;
  }
  
  deleteDofArrays();

  vecU.clear();
  vecDU.clear();
  vecDDU.clear();
  vecR.clear();
  vecF.clear();
}


void FEStorage::listFEComponents() {
  for (size_t i = 0; i < feComponents.size(); i++) {
    LOG(INFO) << *feComponents[i];
  }
}


void FEStorage::initDofs() {
  TIMED_SCOPE(t, "initDofs");
	if (nNodes() < 1 && nElements() < 1) {
		LOG(FATAL) << "No any nodes or elements";
	}
  
  nodeDofs.initDofTable(nNodes());
  elementDofs.initDofTable(nElements());

	for (uint32 el = 0; el < nElements(); el++) {
    elements[el]->pre();
  }

  for (size_t i = 0; i < mpcCollections.size(); i++) {
    mpcCollections[i]->pre();
    mpcCollections[i]->registerMpcsInStorage();
  }

  // Total number of dofs (only registered by elements)
	_nDofs = elementDofs.getNumberOfUsedDofs() + nodeDofs.getNumberOfUsedDofs();
  CHECK(_nDofs);

  auto udofs = getUniqueNodeDofTypes();
  if (udofs.size()) {
    std::stringstream ss;
    ss << "Types of nodal DoFs:";
    for (auto& type : udofs)
      ss << " " << Dof::dofType2label(type);
    LOG(INFO) << ss.str();
  }
  LOG(INFO) << "Number of nodal DoFs: " << nodeDofs.getNumberOfUsedDofs();

  udofs = getUniqueElementDofTypes();
  if (udofs.size()) {
    std::stringstream ss;
    ss << "Types of element DoFs:";
    for (auto& type : udofs)
      ss << " " << Dof::dofType2label(type);
    LOG(INFO) << ss.str();
  }
  LOG(INFO) << "Number of element DoFs: " << elementDofs.getNumberOfUsedDofs();
}


void FEStorage::assignEquationNumbers() {
  // _nUnknownDofs - number of Dof need to be found on every step
	_nUnknownDofs = nDofs() - nConstrainedDofs();
  CHECK(_nUnknownDofs);

	uint32 next_eq_solve = nConstrainedDofs() + 1;
	uint32 next_eq_const = 1;

  // TODO: we can try to do numbering in more efficient way if will loop over element nodes  (as it
  // does in buildK() procedure)
  for (uint32 i = 1; i <= nElements(); i++) {
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

  for (uint32 i = 1; i <= nNodes(); i++) {
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

  assert(next_eq_const - 1 == nConstrainedDofs());
  assert(next_eq_solve - 1 == nDofs());

  for (auto& mpc : mpcs) {
    assert(mpc->eq.size());
    mpc->eqNum = next_eq_solve++;
  }

  assert(next_eq_solve - 1 == nDofs() + nMpc());

	LOG(INFO) << "DoFs = " << nDofs() << ", constrained DoFs = " <<  nConstrainedDofs() << ", MPC eq. = "
      <<  nMpc() << ", TOTAL eq. = " << nUnknownDofs() + nMpc();
}

void FEStorage::initSolutionData () {
  TIMED_SCOPE(t, "initSolutionData");
  
  // We need to know topology of the mesh in order to determine SparsityInfo for sparse matrices
  learnTopology();


  // In nla3d solution procedure there are 3 distinguish types of unknowns. First one "c" -
  // constrained DoFs, and consequently known at solution time. Second one "s" - unknown DoFs (need
  // to be solved).  And last one "lambda" - Lagrange multipliers for applied MPC constraints.
  // Equations numbers assigned in a such manner that constrained DoFs equations go first, thus the
  // global System of Linera Algebraic Equations looks llike this:
  //
  //  |  Kcc  |    KcsMPCc   |   |  Uc   |   | Fc |   | Rc |
  //  |-------|--------------|   |-------|   |    |   |    |
  //  |       |              | * |  Us   | = | Fs | + | Rs |
  //  |KcsMPCc|    KssMPCs   |   |-------|   |    |   |    |
  //  |  ^T   |              |   |  Ul   |   | Fl |   | Rl |
  //  
  //  1. But really only
  //
  //  |              |   |  Us   |   |    |
  //  |    KssMPCs   | * |-------| = |RHS |
  //  |              |   |  Ul   |   |    |
  //
  //  will be be solved by eq. solver. 
  //
  //  where:
  //
  //  |    |   | Fs |   | Rs |    |         |   |    |
  //  |RHS | = |    | + |    | -  |KcsMPCc^T| * | Uc |
  //  |    |   | Fl |   | Rl |    |         |   |    |
  //
  //  2. Then to restore reaction loads:
  //
  //  |    |   |    |   |     |   |    |   |       |   | Us |
  //  | Rc | =-| Fc | + | Kcc | * | Uc | + |KcsMPCc| * |    |
  //  |    |   |    |   |     |   |    |   |       |   | Ul |

  matK = new BlockSparseSymMatrix<2>({nConstrainedDofs(), nUnknownDofs() + nMpc()});

  if (transient) {
    // share sparsity info with K matrices
    matC = new BlockSparseSymMatrix<2>(matK);

    matM = new BlockSparseSymMatrix<2>(matK);
  }


	// initialize solutions vectors
  vecU.reinit(nDofs()+nMpc());
  if (transient) {
    vecDU.reinit(nDofs()+nMpc());
    vecDDU.reinit(nDofs()+nMpc());
  }
  vecR.reinit(nDofs()+nMpc());
  vecF.reinit(nDofs()+nMpc());

	if (!material) {
		LOG(WARNING) << "FEStorage::initializeSolutionData: material isn't defined";
	}


  // Need to restore non-zero entries in Sparse Matrices based on mesh topology and registered Dofs
  // As far as we know from topology which elements are neighbors to each other we can estimate
  // quantity and positions of non-zero coef. in Sparse Matrix matK. Other matrices matC, matK will
  // have the same sparsity as stiffness matrix matK.

  for (uint32 nn = 1; nn <= nNodes(); nn++) {
    auto nn_dofs = nodeDofs.getEntityDofs(nn);

    for (auto& en : topology[nn-1]) {
      // register element dofs to node nn
      auto en_dofs = elementDofs.getEntityDofs(en);
      for (auto d1 = en_dofs.first; d1 != en_dofs.second; d1++)
        for (auto d2 = nn_dofs.first; d2 != nn_dofs.second; d2++)
          addEntryK(d1->eqNumber, d2->eqNumber);

      // cycle over element en Nodes and register nn vs nn2 nodes dofs
      for (uint16 enn = 0; enn < getElement(en).getNNodes(); enn++) {
        uint32 nn2 = getElement(en).getNodeNumber(enn);
        if (nn2 < nn) continue;
        // register node nn2 dofs to node nn dofs
        auto nn2_dofs = nodeDofs.getEntityDofs(nn2);
        for (auto d1 = nn2_dofs.first; d1 != nn2_dofs.second; d1++)
          for (auto d2 = nn_dofs.first; d2 != nn_dofs.second; d2++)
            addEntryK(d1->eqNumber, d2->eqNumber);
      }
    }
  }

  // register element dofs vs element dofs
  for (uint32 en = 1; en <= nElements(); en++) {
    auto en_dofs = elementDofs.getEntityDofs(en);
    for (auto d1 = en_dofs.first; d1 != en_dofs.second; d1++)
      for (auto d2 = en_dofs.first; d2 != en_dofs.second; d2++)
        addEntryK(d1->eqNumber, d2->eqNumber);
  }

  // register MPC coefficients
  for (auto& mpc : mpcs) {
    assert(mpc->eq.size());
    uint32 eq_num = mpc->eqNum;
    for (auto& term : mpc->eq) {
    uint32 eq_j = getNodeDofEqNumber(term.node, term.node_dof);
      addEntryMPC(mpc->eqNum, eq_j);
    }
  }

  // compress sparsity info. After that we can't add new position in sparse matrices
  matK->compress();

  if (transient) {
    matC->compress();
    matM->compress();
  }
}


void FEStorage::printDofInfo(std::ostream& out) {
  for (uint32 en = 1; en <= nElements(); en++) {
    auto en_dofs = elementDofs.getEntityDofs(en);
    for (auto d1 = en_dofs.first; d1 != en_dofs.second; d1++)
      out << "E" << en << ":" << Dof::dofType2label(d1->type) << " eq = " << d1->eqNumber 
           << " constrained = " << d1->isConstrained << endl;
  }

  for (uint32 nn = 1; nn <= nNodes(); nn++) {
    auto nn_dofs = nodeDofs.getEntityDofs(nn);
    for (auto d1 = nn_dofs.first; d1 != nn_dofs.second; d1++)
      out << "N" << nn << ":" << Dof::dofType2label(d1->type) << " eq = " << d1->eqNumber 
           << " constrained = " << d1->isConstrained << endl;
  }

}


void FEStorage::updateResults() {
  TIMED_SCOPE(t, "updateSolutionResults");
  // calculate element's update procedures (calculate stresses, strains, ..)
  for (uint32 el = 0; el < nElements(); el++) {
    elements[el]->update();
  }
}


void FEStorage::learnTopology() {
  topology.clear();
  topology.assign(nNodes(), std::set<uint32>());
  for (uint32 en = 1; en <= nElements(); en++) {
    for (uint16 nn = 0; nn < getElement(en).getNNodes(); nn++) {
      uint32 noden = getElement(en).getNodeNumber(nn);
      topology[noden-1].insert(getElement(en).getElNum());
    }
  }
}


} // namespace nla3d 

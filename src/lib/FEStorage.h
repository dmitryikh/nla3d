// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include <list>
#include "sys.h"
#include "materials/MaterialFactory.h"
#include "elements/ElementFactory.h"
#include "math/BlockSparseMatrix.h"
#include "FEComponent.h"
#include "Mpc.h"
 
namespace nla3d {

class Element;
class Node;
class Dof;
class ElementFactory;


// FEStorage - heart of the nla3d program, it contains all FE data: 
//   * nodes array `nodes`, elements array `elements`,
//   * nodal DoFs array `elementDofs`, element DoFs array `nodeDofs`,
//   * solution matrices and vectors (`matK`, `vecU`, `vecF`, `vecR` and others),
//   * mesh topology (connectivity) info `topology`,
//   * mpc equations `mpcs`.
//
// FEStorage provides methods to fill and work with FE data described above:
//   * methods to add nodes and elements to FEStorage
//   * methods to initialize solution data and provide it to FESolver
//   * methods to assemble global system of equations in form of: M * DDU + C * DU + K * U = F + R
//   * a lot of getters to get FE data
//
// NOTE: Numbering. Node, element, MPC equation, Equation numbers are started from 1.
// NOTE: FEStorage works with dense node, element arrays. That means that element/node numbers
// should start from 1 and up to nElements()/nNodes().
class FEStorage {
public:
	FEStorage();
	~FEStorage();

  // A pointer to Material instance. nla3d supports only one material for FE model.  External code
  // should create material instance and pass it to FEStorage. Then FEStorage takes a control on
  // material. The material will be deleted in ~FEStorage().
  // NOTE: material conception should be implement in Element class. In later updates
  // FEStorage::Material will be deleted.
  Material* material = nullptr;

  // *** Operations for preparing the global system of equations ***
  //
  // These functions add a value to an entry of global system of equations matrices (stiffness,
  // damping, inertia).
  // For example: addValueK(32, Dof::UX, 42, Dof::UZ, 0.32) will add 0.32 to a corresponding
  // entry in a global stiffness matrix.
  // NOTE: node > 0
  // NOTE: use addValueK(uint32 nodei, ...) to add nodal Dof vs nodal Dof value
  // or use addValueK(uint32 eqi, ..) to add value for equation number i vs equation number j (more general case)
  void addValueK(uint32 nodei, Dof::dofType dofi, uint32 nodej, Dof::dofType dofj, double value);
  void addValueK(uint32 eqi, uint32 eqj, double value);

  // for damping matrix
  void addValueC(uint32 nodei, Dof::dofType dofi, uint32 nodej, Dof::dofType dofj, double value);
  void addValueC(uint32 eqi, uint32 eqj, double value);

  // for mass matrix
  void addValueM(uint32 nodei, Dof::dofType dofi, uint32 nodej, Dof::dofType dofj, double value);
  void addValueM(uint32 eqi, uint32 eqj, double value);

  // Add value to a matrix of MPC coefficients for DoFs.
  // NOTE: MPC equations can work only with vecU values, excluding derivatives values
  // vecDU, vecDDU.
  void addValueMPC(uint32 eq_num, uint32 nodej, Dof::dofType dofj, double coef);
  void addValueMPC(uint32 eq_num, uint32 eqj, double coef);

  // Add value to vecF (term of RHS of global equations system representing internal element loads)
  // NOTE: use addValueF(uint32 nodei, ...) to add value to vecF for nodal Dof
  // or use addValueF(uint32 eqi, ..) to add value for equation number i (more general case)
  // NOTE: common way is that element's assemble procedure contributes to vecF by mean of these
  // functions
  void addValueF(uint32 nodei, Dof::dofType dofi, double value);
  void addValueF(uint32 eqi, double value);

  // Add value to vecR (term of RHS of global equations system representing external loads)
  // NOTE: common way is that FESolver or user code contribute to vecR by mean of these
  // functions
  void addValueR(uint32 nodei, Dof::dofType dofi, double value);
  void addValueR(uint32 eqi, double value);

  // fill with zeros the matrices of global system of equations
	void zeroK();
	void zeroC();
	void zeroM();
  // fill with zeros the vector vecF
	void zeroF();

  // get the K, C, M matrices (stiffness, damping, inertia) of the global system of equations.  In current realization the
  // matrix is symmetric, but semi-positive defined, because of MPC equations (zeros on diagonals in
  // Lagrangian columns x rows)
  math::BlockSparseSymMatrix<2>* getK();
  math::BlockSparseSymMatrix<2>* getC();
  math::BlockSparseSymMatrix<2>* getM();

  // get RHS vectors of the global system of equations
  math::dVec* getF();
  math::dVec* getR();

  // get DoF values vectors
  math::dVec* getU();
  math::dVec* getDU();
  math::dVec* getDDU();

  // Procedure of filling global equations system matrices and RHS vectors with actual values.
  // if isTransient() == bool then matC, matM are also assembled.
  void assembleGlobalEqMatrices();

  // getters to get numbers of different entities stored in FEStorage
	uint32 nNodes();
	uint32 nElements();
	uint32 nDofs();
	uint32 nUnknownDofs();
	uint32 nConstrainedDofs();
	uint32 nMpc();

  // if isTransient() == bool then FEStorage initialize matC, matM, vecDU, vecDDU along with matK,
  // vecU. 
  void setTransient(bool _transient);
  bool isTransient();

  // Operations with DoFs
  //
  // Registation of DoFs is a key moment in nla3d. Every element (and other entities like MPC
  // collections) should register its DoFs which the element will use. 
  // For example: addNodeDof(32, Dof::UX) means that DoF UX for node 32 will be used in a global system of
  // equations. Unregistered DoFs will not participate in global equations system. This mechanism
  // leads to ability to have different DoFs in different nodes. 
  // NOTE: `node` > 0
  void addNodeDof(uint32 node, std::initializer_list<Dof::dofType> _dofs);
  // NOTE: `element` > 0
  void addElementDof(uint32 el, std::initializer_list<Dof::dofType> _dofs);
  // functions to determine how many uniques DoF types are used for nodes and elements
  std::set<Dof::dofType> getUniqueNodeDofTypes();
  std::set<Dof::dofType> getUniqueElementDofTypes();
  // return true if the corresponding DoF was added in nodeDofs or elementDofs.
  bool isElementDofUsed(uint32 el, Dof::dofType dof);
  bool isNodeDofUsed(uint32 node, Dof::dofType dof);
  // get the number of equation in the global eq. system for particular DoF
	uint32 getElementDofEqNumber(uint32 el, Dof::dofType dof);
	uint32 getNodeDofEqNumber(uint32 node, Dof::dofType dof);
  // get a pointer to Dof class for particular DoF
  Dof* getElementDof(uint32 el, Dof::dofType dof);
  Dof* getNodeDof(uint32 node, Dof::dofType dof);
  // get the obtained value for a particular DoF. The value is a solution to the present moment
  double getNodeDofSolution(uint32 node, Dof::dofType dof);
  double getElementDofSolution(uint32 el, Dof::dofType dof);
  // as far as constrained (fixed) DoFs have special treatment in FEStorage outer code (commonly
  // FESolver) should register all constrained DoFs by mean of these functions (before calling
  // assignEquationNumbers()) 
  void setConstrainedNodeDof(uint32 node, Dof::dofType dtype);
  void setConstrainedElementDof(uint32 el, Dof::dofType dtype);
  // Get value of reaction for particular DoFs. Reactions are defined only for constrained DoFs.  If
  // a reaction value is asked for not constrained DoFs, a value of external concentrated load is
  // returned.
  // NOTE: actually this functions just return vecR[eq-1] value, for eq  <= nConstrainedDofs() it's
  // reaction load, but for eq > nConstrainedDofs() it's external DoF's concentrated load.
	double getReaction(uint32 node, Dof::dofType dof);
	double getReaction(uint32 eq);

  // get an instance of Material class
	Material* getMaterial();

  // get an instance of particular node
  // NOTE: `_nn` > 0
	Node& getNode(uint32 _nn);
  // function fills node_ptr with pointers to Node classes for element el.
  // Calling side should reserve a space for an array of pointers node_ptr.
  // size of node_ptr should be at least Element::getNNodes()
  // NOTE: `el` > 0
	void getElementNodes(uint32 el, Node** node_ptr);
  // The function return a spatial position of node `n` in either initial state (`deformed` = false)
  // or deformed state (`deformed` = true);
  // A calling side should take care about memory allocation for ptr. Size of ptr should be at least
  // 3 as the function always return 3-dimensional coordinates.
	void getNodePosition(uint32 n, double* ptr, bool deformed = false);
  // NOTE: `_en` > 0
	Element& getElement(uint32 _en);
  template<typename ET>
  ET& getElement(uint32 _en);
  // get a FEComponent instance by registration number
  // NOTE: `i` > -1
  FEComponent* getFEComponent(size_t i);
  // get a FEComponent instance by component name
  FEComponent* getFEComponent(const std::string& name);

  // storing operations
  //
  // Add `mpc` to `mpcs` array. FEStorage will delete Mpc instances by itslef.
	void addMpc(Mpc* mpc);
  // Add `mpcCol` to `mpcCollections` array. FEStorage will delete MpcCollection instances by itslef.
	void addMpcCollection(MpcCollection* mpcCol);
  // Add `comp` to `feComponents` array. FEStorage will delete FEComponent instances by itslef.
  void addFEComponent(FEComponent* comp);
  // The function add a node to the nodes table.
  // TODO: this is inefficient functions because we can't ensure memory localization for all nodes
  // instanses in `nodes` array
  void addNode(Node* node);
  // The function creates an vector `nodes` with dynamically allocated Node instances inside.
  // Numbers of newly created elements pass back to the caller.  
  // NOTE: nodes will be created with default Node() constructor (node coordinates 0, 0, 0). 
  // NOTE: new nodes are concantenated to old nodes which already were in FEStorage 
	std::vector<uint32> createNodes(uint32 nn); 
  // add an element to the element array `elements`.
  // TODO: this is inefficient functions because we can't ensure memory localization for all
  // elements instances in `elements` array
  void addElement(Element* el);
  // The function creates a vector with dynamically allocated Element instances inside.  Particular
  // realization of abstract Element class is chosen from elType variable by mean of ElementFactory
  // class (see elements/ElemenetFactory.h). Numbers of newly created elements pass back to the
  // caller. 
  // NOTE: elements will be created with default constructor. User code should then fill node
  // numbers (by Element::getNodeNumber(..)) and other stuff related to particular realization of
  // Element class.
  std::vector<uint32> createElements(uint32 en, ElementType elType);
  // template version of function described above. User code should provide example Element entity
  // of particular class. 
  template<typename T>
  std::vector<uint32> createElements(uint32 _en, T example);

  // delete procedures
  //
  // delete all FE model things: elements, nodes, MPCs, MPC Collections, FE components, topology
  // info.
	void deleteMesh();
  // delete nodes table: delete all dynamically allocated Node instances,
  // and clear vector of pointers `nodes`.
  void deleteNodes();
  // delete elements table: delete all dynamically allocated Element instances,
  // and clear vector of pointers `elements`.
  void deleteElements();
  // delete MPC list: delete all dynamically allocated Mpc instances,
  // and clear vector of pointers `mpcs`.
  void deleteMpcs();
  void deleteMpcCollections();
  void deleteFeComponents();
  // delete elementDofs and nodeDofs arrays (and free the memory)
  void deleteDofArrays();
  // delete all data allocated in FEStorage::initSolutionData(), delete DoFs arrays
  void deleteSolutionData();

  // print procedures
  void listFEComponents();
  // for debug purpose only. Be carefully, this is output intensive..
  void printDofInfo(std::ostream& out);

  // solution procedures
  // These functions perform all preparing steps before FESolver will solve the problem:
  // After this procedure one can't add more elements, nodes, boundary conditions, MPCs, ..
  // 1. Fill DoFs tables by calling Element->pre() and MpcCollection->pre() to register needed DoFs.
	void initDofs();
  // 2. Outer code should call FEStorage::setConstrained[Node/Element]Dof(..) to set the particular
  // DoFs to be constrained (fixed).
  // 3. Assign number for globals system equations by calling assignEquationNumbers();
  void assignEquationNumbers();
  // 4. Allocate all solution data structures: matK/C/M, vecU/DU/DDU/F/R
	void initSolutionData();

  // After global equations system is solved and vecU/DU/DDU/R is updated with appropriate values
  // FESolver should call this procedure to update element solution data
  // NOTE: actually Element::update() is called
	void updateResults();

private:
  // fill `topology` data based on the current mesh (Element::nodes numbers)
  void learnTopology();

  // these functions are used to train sparsity info for matK/C/M
  // provide info that entry (eqi, eqj) are not zero
  // should be called before matK->compressed()
  void addEntryK(uint32 eqi, uint32 eqj);
  void addEntryMPC(uint32 eq_num, uint32 eqj);

  // Total number of DoFs (registered by FEStorage::add[Node/Element]Dof(..))
	uint32 _nDofs = 0;
  // Number of constrained (fixed) DoFs values (registered by
  // FEStorage::setConstrained[Node/Element]Dof())
	uint32 _nConstrainedDofs = 0;
  // _nUnknownDofs - number of Dof need to be found on every step (excluding Lagrangian lambdas)
	uint32 _nUnknownDofs;

  // Array of elements. Elements are created by the FEStorage::createElements(..) function.  It's
  // important that elements in nla3d have consecutive numbering. That means that they have numbers
  // in [1; nElements()].
  // NOTE: Elements are deleted by deleteElements() function.
	std::vector<Element*> elements;

  // Array of nodes. Nodes are created by FEStorage::createNodes(..). Nodes are created with default
  // coordinates (0,0,0). Nodes have consecutive numbering. They have numbers in [1; nNodes()]. 
	std::vector<Node*> nodes;

  // List of MPC equations. List is populated by FEStorage::addMpc(..). An instance of Mpc class is
  // created outside of FEStorage class. But after addMpc(..) function FEStorage takes control on
  // the instance of Mpc. MPCs can be deleted by FEStorage::deleteMpcs(). Of course, they are
  // deleted in ~FEStorage() destructor. 
	std::list<Mpc*> mpcs;

  // FE components - lists consist of numbers of entities. Here is the support of nodal components
  // and element components. This components can be used to apply BCs and/or MPCs. FE component is
  // attached to FE Storage by FEStorage::addFEComponent(..) function. All stored components can be
  // listed by FEStorage::listFEComponents(..). One can access to the particular component by
  // methods FEStorage:: getFEComponent(..). 
  std::vector<FEComponent*> feComponents;

  // An array of MpcCollection. Some number of MPCs can be grouped together because of they were
  // produced by a single rule. For example, in FE model can rigid body MPC equations between the
  // master node and slave nodes.  The relation between the master node and slave nodes described by
  // a group of MPC equations. All of them are grouped in MpcCollection. See Mpc.h for more details.
  // Here is a rule: Mpc defined inside of MpcCollection should be added to FEStorage::mpcs by
  // addMpc(..).  That means that Mpc instances are dynamically created in MpcCollection, but then
  // FEStorage class takes control on it.
  std::vector<MpcCollection*> mpcCollections;

  // Stiffness matrix. It's BlockSparseSymMatrix because of block(1) and block(1, 2) related to
  // equations for constrained (fixed) DoFs and will not be used in equation solver. Instead of this
  // only block(2) should be solved by eq. solver.
  // NOTE: Such separation is possible because FEStorage gives to constrained DoFs equation numbers
  // from 1 to nConstrainedDofs(). Unknown DoFs has numbers from nConstrainedDofs() + 1 to
  // nConstrainedDofs() + nUnknownDofs(). 
  math::BlockSparseSymMatrix<2>* matK = nullptr;

  // matrices for transient analysis:
  // NOTE: matK, matC, matM have the same SparsityInfo for fast math like this : matKmod = c1 * matK + c2 *
  // matC + c3 * matM.
  // Damping matrix
  math::BlockSparseSymMatrix<2>* matC = nullptr;
  // Inertia matrix
  math::BlockSparseSymMatrix<2>* matM = nullptr;

  // A found values for DoFs. vecU has size of nDofs + nMpc(). FESolver is responsible of filling
  // this vector.
  // The whole vector could be divided into different parts: values for constrained DoFs in range
  // [0; nConstrainedDofs()), values for unknown DoFs [nConstrainedDofs(); nDofs()), values for
  // Mpc's Lagrangian lambdas [nDofs(); nDofs() + nMpc()). 
  math::dVec vecU;

  // Found first time derivatives for DoF's values. FESovler is responsible of filling this vector.
  math::dVec vecDU;

  // Found second time derivatives for DoF's values. FESovler is responsible of filling this vector.
  math::dVec vecDDU;

  // vector stores external concentrated loads. vecR has size nDofs + nMpc(). FESolver is responsible of filling
  // this vector.
  // The whole vector could be divided into different parts: values of reactions for constrained DoFs in range
  // [0; nConstrainedDofs()), values of external loads for unknown DoFs [nConstrainedDofs(); nDofs()), for
  // Mpc's Lagrangian lambdas [nDofs(); nDofs() + nMpc()) values should be alway remain zeros. 
  math::dVec vecR;

  // internal element loads. vecF has size nDofs + nMpc(). FEStorage::assembleGlobalEqMatrices() is
  // responsible of filling this vector.
  math::dVec vecF;

  // * In nla3d every DoF is represented by `Dof` class.
  // * DoF can be nodal or for an element.
  // * `elementDofs` and `nodeDofs` are DofCollection 
  //    DofCollection keeps Dof objects and give it by (node, dof) pair
  // * Elements (or another entities) should register DoFs by using
  //    FEStorage::add[Node/Element]Dof(node, dof)
  DofCollection elementDofs;
  DofCollection nodeDofs;

  // topology array. It stores a set of element numbers [el1, el2, el3..] attached to the node `n`:
  // topology[n-1] = [el1, el2, el3..]
  std::vector<std::set<uint32> > topology;

  // if transient is true that means that assembleGlobalEqMatrices() should assemble M and C
  // matrices too
  bool transient = false;
};


inline void FEStorage::addValueK(uint32 nodei, Dof::dofType dofi, uint32 nodej, Dof::dofType dofj, double value) {
	uint32 rowEq = getNodeDofEqNumber(nodei, dofi);
	uint32 colEq = getNodeDofEqNumber(nodej, dofj);
  addValueK(rowEq, colEq, value);
}


inline void FEStorage::addValueK(uint32 eqi, uint32 eqj, double value) {
  // eqi - row equation 
  // eqj - column equation
  matK->addValue(eqi, eqj, value);
}


inline void FEStorage::addValueC(uint32 nodei, Dof::dofType dofi, uint32 nodej, Dof::dofType dofj, double value) {
	uint32 rowEq = getNodeDofEqNumber(nodei, dofi);
	uint32 colEq = getNodeDofEqNumber(nodej, dofj);
  addValueC(rowEq, colEq, value);
}


inline void FEStorage::addValueC(uint32 eqi, uint32 eqj, double value) {
  // eqi - row equation 
  // eqj - column equation
  matC->addValue(eqi, eqj, value);
}


inline void FEStorage::addValueM(uint32 nodei, Dof::dofType dofi, uint32 nodej, Dof::dofType dofj, double value) {
	uint32 rowEq = getNodeDofEqNumber(nodei, dofi);
	uint32 colEq = getNodeDofEqNumber(nodej, dofj);
  addValueM(rowEq, colEq, value);
}


inline void FEStorage::addValueM(uint32 eqi, uint32 eqj, double value) {
  // eqi - row equation 
  // eqj - column equation
  matM->addValue(eqi, eqj, value);
}


inline void FEStorage::addValueMPC(uint32 eq_num, uint32 nodej, Dof::dofType dofj, double coef) {
	uint32 colEq = getNodeDofEqNumber(nodej, dofj);
  addValueMPC(eq_num, colEq, coef);
}


inline void FEStorage::addValueMPC(uint32 eq_num, uint32 eqj, double coef) {
  matK->addValue(eq_num, eqj, coef);
}


inline void FEStorage::addValueF(uint32 nodei, Dof::dofType dofi, double value) {
	uint32 rowEq = getNodeDofEqNumber(nodei, dofi);
  addValueF(rowEq, value);
}


inline void FEStorage::addValueF(uint32 eqi, double value) {
  assert(eqi > 0);
	assert(eqi <= vecF.size());
	assert(eqi <= nDofs() + nMpc());
	vecF[eqi - 1] += value;
}


inline void FEStorage::addValueR(uint32 nodei, Dof::dofType dofi, double value) {
	uint32 rowEq = getNodeDofEqNumber(nodei, dofi);
  addValueR(rowEq, value);
}


inline void FEStorage::addValueR(uint32 eqi, double value) {
  assert(eqi > 0);
	assert(eqi <= vecR.size());
	assert(eqi <= nDofs() + nMpc());
	vecR[eqi - 1] += value;
}


inline void FEStorage::zeroK() {
	assert(matK);
  matK->zero();
}


inline void FEStorage::zeroC() {
	assert(matC);
  matC->zero();
}


inline void FEStorage::zeroM() {
	assert(matM);
  matM->zero();
}


inline void FEStorage::zeroF() {
  vecF.zero();
}


inline math::BlockSparseSymMatrix<2>* FEStorage::getK() {
	assert(matK->isCompressed());
	return matK;
}


inline math::BlockSparseSymMatrix<2>* FEStorage::getC() {
	assert(matC->isCompressed());
	return matC;
}


inline math::BlockSparseSymMatrix<2>* FEStorage::getM() {
	assert(matM->isCompressed());
	return matM;
}


inline math::dVec* FEStorage::getF() {
  assert(vecF.size() > 0);
  return &vecF;
}


inline math::dVec* FEStorage::getU() {
  assert(vecU.size() > 0);
  return &vecU;
}


inline math::dVec* FEStorage::getDU() {
  assert(vecDU.size() > 0);
  return &vecDU;
}


inline math::dVec* FEStorage::getDDU() {
  assert(vecDDU.size() > 0);
  return &vecDDU;
}


inline math::dVec* FEStorage::getR() {
  assert(vecR.size() > 0);
  return &vecR;
}

inline uint32 FEStorage::nDofs() {
  return _nDofs;
}


inline uint32 FEStorage::nConstrainedDofs() {
  return _nConstrainedDofs;
}

inline uint32 FEStorage::nUnknownDofs() {
  return _nUnknownDofs;
}


inline uint32 FEStorage::nMpc() {
  return static_cast<uint32> (mpcs.size());
}


inline uint32 FEStorage::nNodes () {
  return static_cast<uint32> (nodes.size());
}

inline uint32 FEStorage::nElements () {
  return static_cast<uint32> (elements.size());
}


inline void FEStorage::setTransient(bool _transient) {
  transient = _transient;
}


inline bool FEStorage::isTransient() {
  return transient;
}

inline void FEStorage::addNodeDof(uint32 node, std::initializer_list<Dof::dofType> _dofs) {
  assert(nodeDofs.getNumberOfEntities() > 0);
  nodeDofs.addDof(node, _dofs);
}

inline void FEStorage::addElementDof(uint32 el, std::initializer_list<Dof::dofType> _dofs) {
  assert(elementDofs.getNumberOfEntities() > 0);
  elementDofs.addDof(el, _dofs);
}


inline std::set<Dof::dofType> FEStorage::getUniqueNodeDofTypes() {
  return nodeDofs.getUniqueDofTypes();
}


inline std::set<Dof::dofType> FEStorage::getUniqueElementDofTypes() {
  return elementDofs.getUniqueDofTypes();
}


inline bool FEStorage::isElementDofUsed (uint32 el, Dof::dofType dof) {
  assert(elementDofs.getNumberOfEntities() > 0);
  return elementDofs.isDofUsed(el, dof);
}


inline bool FEStorage::isNodeDofUsed (uint32 node, Dof::dofType dof) {
  assert(nodeDofs.getNumberOfEntities() > 0);
  return nodeDofs.isDofUsed(node, dof);
}


inline uint32 FEStorage::getElementDofEqNumber(uint32 el, Dof::dofType dof) {
  return getElementDof(el, dof)->eqNumber;
}


inline uint32 FEStorage::getNodeDofEqNumber(uint32 node, Dof::dofType dof) {
  return getNodeDof(node, dof)->eqNumber;
}


inline Dof* FEStorage::getElementDof (uint32 el, Dof::dofType dof) {
  assert(elementDofs.getNumberOfUsedDofs() > 0);
  return elementDofs.getDof(el, dof);
}

inline Dof* FEStorage::getNodeDof (uint32 node, Dof::dofType dof) {
  assert(nodeDofs.getNumberOfUsedDofs() > 0);
  return nodeDofs.getDof(node, dof);
}


inline double FEStorage::getNodeDofSolution (uint32 node, Dof::dofType dof) {
  assert(vecU.size() > 0);
  return vecU[getNodeDofEqNumber(node, dof) - 1];
}


inline double FEStorage::getElementDofSolution (uint32 el, Dof::dofType dof) {
  assert(vecU.size() > 0);
  return vecU[getElementDofEqNumber(el, dof) - 1];
}


inline Node& FEStorage::getNode(uint32 _nn) {
	assert(_nn> 0 && _nn <= nNodes());
	return *nodes[_nn-1];
}


inline Element& FEStorage::getElement(uint32 _en) {
	assert(_en <= nElements());
	return *(elements[_en-1]);
}


template<typename ET>
inline ET& FEStorage::getElement(uint32 _en) {
  return dynamic_cast<ET&>(getElement(_en));
}


inline void FEStorage::addEntryK(uint32 eqi, uint32 eqj) {
  // eqi - row equation 
  // eqj - column equation
  matK->addEntry(eqi, eqj);
}


inline void FEStorage::addEntryMPC(uint32 eq_num, uint32 eqj) {
	assert(matK);
	assert(eq_num > nConstrainedDofs() + nUnknownDofs());
  assert(eq_num <= nConstrainedDofs() + nUnknownDofs() + nMpc());
  assert(eqj > 0 && eqj <= nConstrainedDofs() + nUnknownDofs());
  matK->addEntry(eq_num, eqj);
}


} // namespace nla3d 

// 'dirty' hack to avoid include loops (element-vs-festorage)
#include "elements/element.h"

namespace nla3d {
template<typename T>
std::vector<uint32> FEStorage::createElements(uint32 _en, T example) {
  //TODO: catch if not enough memory
  std::vector<uint32> newIndexes;
  newIndexes.reserve(_en);
  uint32 nextNumber = elements.size() + 1;

  elements.reserve(elements.size() + _en);

  for (uint32 i = 0; i < _en; i++) {
    T* els = new T;
    elements.push_back(els);
  }

  for (uint32 i = nextNumber; i <= elements.size(); i++) {
    //access elNum protected values as friend
    elements[i - 1]->elNum = i;
    elements[i - 1]->storage = this;
    newIndexes.push_back(i);
  }
  return newIndexes;
}



} // namespace nla3d 

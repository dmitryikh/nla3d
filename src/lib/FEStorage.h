// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include <list>
#include "sys.h"
#include "materials/MaterialFactory.h"
#include "elements/ElementFactory.h"
#include "math/SparseMatrix.h"
#include "FEComponent.h"
#include "Mpc.h"
 
namespace nla3d {

//pre-defines
class Element;
class Node;
class Dof;
class ElementFactory;


// class FEStorage - heart of the NLA programm, it contains all data to keep in memory: 
//  nodes, elements, DoFs, internal data to map equation numbers and DoFs, 
//	solution staff (global matrixes and vectors), processors.

class FEStorage {
public:
	FEStorage();
	~FEStorage();

  // A pointer to Material instance. nla3d supports only one material for FE model.  Thats why this
  // here is just a pointer to a material. External code should create material instance and pass it
  // to FEStorage. Then FEStorage takes a control on material. The material will be deleted in
  // ~FEStorage().
  Material* material;

  // Operations for preparing the global system of equations
  //
  // These functions add a value to an element of global system of equations matrix (stiffness matrix).
  // For example: addValueK(32, Dof::UX, 42, Dof::UZ, 0.32) will add 0.32 to a corresponding
  // element in a global stiffness matrix.
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
  void addValueMPC(uint32 eq_num, uint32 nodej, Dof::dofType dofj, double coef);
  void addValueMPC(uint32 eq_num, uint32 eqj, double coef);
  // Add value to RHS of global system of equations
  // NOTE: node > 0
  // NOTE: use addValueF(uint32 nodei, ...) to add nodal Dof RHS value
  // or use addValueF(uint32 eqi, ..) to add value for equation number i (more general case)
  void addValueF(uint32 nodei, Dof::dofType dofi, double value);
  void addValueF(uint32 eqi, double value);

  // fill with zeros the matrix of global system of equations. As far as global system of equations consist of
  // different blocks, zeroK() zeros next sparse matrices: KssMPCsT, Kcc, Kcs, MPCc 
	void zeroK();
	void zeroC();
	void zeroM();
  // fill with zeros the vector of RHS (right hand side), solutionDeltaValues vector
  // and mpcConstantValues.
	void zeroF();

  // get the K matrix (stiffness) of the global system of equations.  In current realization the
  // matrix is symmetric, but not positive defined, because of MPC equations (zeros on diagonals in
  // Lagrangian columns x rows) 
  math::SparseSymMatrix* getK();
  math::SparseSymMatrix* getC();
  math::SparseSymMatrix* getM();

  // get right hand side of the global system of equations. 
	double* getF();
  // get vector of Unknowns. Before solution procedure the vector is set to zero.
  // After the solution procedure the vector contain the solution of current equilibrium iteration
	double* getGlobalEqUnknowns();
  // Procedure of filling global equations system with coefficients.
  void assembleGlobalEqMatrices();

	uint32 getNumberOfNodes ();
	uint32 getNumberOfElements ();
	uint32 getNumberOfUnknownDofs ();
	uint32 getNumberOfConstrainedDofs ();
	uint32 getNumberOfMpcEq ();

  void setTransient(bool _transient);
  bool isTransient();

  // Operations with DoFs
  //
  // registation of DoFs is a key moment in nla3d. Every element (and other things) should register those
  // DoFs which the element will use. 
  // For example: addNodeDof(32, Dof::UX) means that DoF UX for node 32 will be used in a global system of
  // equations. Unregistered DoFs will be eliminated automatically. That means that such DoFs contain only zero 
  // rows and columns. That leads to ability to have different DoFs in different nodes/elements. 
  // node number is always > 0
  void addNodeDof(uint32 node, std::initializer_list<Dof::dofType> _dofs);
  // element number is always > 0
  void addElementDof(uint32 el, std::initializer_list<Dof::dofType> _dofs);
  // functions to determine how many uniques DoF types are used for nodes and elements
  std::set<Dof::dofType> getUniqueNodeDofTypes();
  std::set<Dof::dofType> getUniqueElementDofTypes();
  // return true if the corresponding DoF was registered as used in nodeDofs or elementDofs.
  bool isElementDofUsed(uint32 el, Dof::dofType dof);
  bool isNodeDofUsed(uint32 node, Dof::dofType dof);
  // get the number of equation in the global eq. system for particular DoF
	uint32 getElementDofEqNumber(uint32 el, Dof::dofType dof);
	uint32 getNodeDofEqNumber(uint32 node, Dof::dofType dof);
  // get instance of Dof class for particular DoF
  Dof* getElementDof (uint32 el, Dof::dofType dof);
  Dof* getNodeDof (uint32 node, Dof::dofType dof);
  // get the obtained value for a particular DoF. The value is a solution to the present moment
  // (last converged equilibrium iteration).
  double getNodeDofSolution (uint32 node, Dof::dofType dof);
  double getElementDofSolution (uint32 el, Dof::dofType dof);
  // get the last delta value after solving last converged equilibrium iteration
  double getNodeDofSolutionDelta (uint32 node, Dof::dofType dof);
  double getElementDofSolutionDelta (uint32 el, Dof::dofType dof);
  // Get magnitude of reaction for particular DoFs. Reactions are defined only for constrained DoFs.
  // If a reaction value is asked for not constrained DoFs, a zero value is returned.
	double getReaction(uint32 node, Dof::dofType dof);
  // by passing equation number
	double getReaction(uint32 eq);


  // get an instance of Material class
	Material* getMaterial();
  // get an instance of particular node
  // _nn > 0
	Node& getNode(uint32 _nn);
  // function fills node_ptr with pointers to Node classes for element el.
  // Calling side should reserve a space for an array of pointers node_ptr.
  // size of node_ptr should be at least Element::n_nodes()
  // el > 0
	void getElementNodes (uint32 el, Node** node_ptr);
  // the function return a spatial position of node n in either initial state (deformed = false)
  // or deformed state (deformed = true);
  // A calling side should take care about memory allocation for ptr. Size of ptr should be
  // at least 3 as the function always return 3-dimensional coordinates.
	void getNodePosition(uint32 n, double* ptr, bool deformed = false); //def. or initial pos
  // _en > 0
	Element& getElement(uint32 _en);
  // get a FEComponent instance by registration number
  // i > -1
  FEComponent* getFEComponent(size_t i);
  // get a FEComponent instance by component name
  FEComponent* getFEComponent(const std::string& name);


  // storing operations
  //
  
  // add DoF constraint boundary condition 
  // n > 0 - nodal DoF, n < 0 element DoF.
  void addDofFixation (int32 n, Dof::dofType dof, const double value = 0.0);
  // add DoF load (force) boundary condition 
  // n > 0 - nodal DoF, n < 0 element DoF.
  void addDofLoad(int32 n, Dof::dofType dof, const double value = 0.0);
	void addBoundaryCondition (BC_dof_constraint &bc) ;
	void addBoundaryCondition (BC_dof_force &bc) ;
	void addMpc (Mpc* mpc);
	void addMpcCollection (MpcCollection* mpcCol);
  void addFEComponent (FEComponent* comp);
  // The function add a node to the nodes table.
  void addNode(Node* node);
  // The function creates a vector with dynamically allocated Node instances inside.
  // All previously created Nodes will be deleted. That means that by this function one can
  // reallocate node table. 
  // NOTE: nodes will be created with default Node() constructor (node coordinates 0,0,0). 
	void createNodes (uint32 nn); 
  // add an element to the element table.
  void addElement (Element* el);
  // The function creates a vector with dynamically allocated Element instances inside.
  // Particular realisation of abstract Element class is choosed from elType variable by mean of
  // ElementFactory class (see elements/ElemenetFactory.h). All previously created elements will
  // be deleted. That means that by this function one can reallocate element table. 
  // NOTE: elements will be created with default constructor. Outer procedure should then fill node
  // numbers (by Element::getNodeNumber(..)) and other stuff related to particular realisation of
  // Element class.
	void createElements (uint32 en, ElementType elType);

  // delete procedures
  //
  // delete all FE model things: elements, nodes, boundary conditions, MPCs.
	void deleteMesh ();
  // delete nodes table: delete all dynamically allocated Node instances,
  // and clear vector of pointers (nodes).
  void deleteNodes ();
  // delete elements table: delete all dynamically allocated Element instances,
  // and clear vector of pointers (element).
  void deleteElements ();
  // delete MPC list: delete all dynamically allocated Mpc instances,
  // and clear vector of pointers (mpcs).
  void deleteMpcs();
  void deleteMpcCollections();
  void deleteFeComponents();
  // delete elementDofs and nodeDofs arrays (and free memory)
  void deleteDofArrays();
  // delete all data allocated in FEStorage::initializeSolutionData()
  void deleteSolutionData();

  // print procedures
  void listFEComponents();

  // solution procedures
  // The function where all preparing steps before solution are done.
  // elements[i]->pre() and mpcCollections[j]->pre() register all DoFs.
  // Then all memory needed to store solution data are allocated.
  // After this procedure one can't add more elements, nodes, boundary conditions, MPCs, ..
	bool initializeSolutionData();

  // for debug purpose only. Be carefully, this is output intensive..
  void printDofInfo(std::ostream& out);

  // After global equations system is solved this procedure updates solution data:
  // * update dofValues;
  // * find reactions for constrained DoFs;
  // * elements[i]->update()
	void updateSolutionResults();

  // TODO: this kludge as far as updateTimestepResults and updateSolutionResults should be combined
  // into one function
  void updateTimestepResults();
  
  // The function applies boundary conditions
  // The function receives current normalized time (in range [0.0; 1.0]) and normalized time delta
	void applyBoundaryConditions(double time, double timeDelta);

private:
  void learnTopology();
  void addEntryK(uint32 eqi, uint32 eqj);
  void addEntryMPC(uint32 eq_num, uint32 eqj);

	uint32 numberOfNodes;
	uint32 numberOfElements;
  // Total number of dofs (only registered by FEStorage::registerNodeDof(..) or FEStorage::registerElementDof(..))
	uint32 numberOfDofs;

  // Number of constrained (fixed by BCs) dofs.
	uint32 numberOfConstrainedDofs;
  // numberOfUnknownDofs - number of Dof need to be found on every step
	uint32 numberOfUnknownDofs;
  // Number of MPC equations
	uint32 numberOfMpcEq;

  // Array of elements. Elements are created by the FEStorage::createElements(..) function.
  // It's important that elements in nla3d have consecutive numbering. That means that they have
  // numbers in [1; getNumberOfElements()]. Particular element type is choosen based on 
  // FEStorage::elType variable. NOTE that in FE model cna be used only one type of element.
  // Elements deleted by deleteElements() function.
	std::vector<Element*> elements;

  // Array of nodes. Nodes are created by FEStorage::createNodes(..). Nodes are created with
  // default coordinates (0,0,0). Nodes have consecutive numbering. They have numbers
  // in [1; getNumberOfNodes()]. 
	std::vector<Node*> nodes;

  // List of force boundary conditions. This is old fashion realisation and will be changed soon.
  // List is populated by FEStorage::addBoundaryCondition(..).
	std::list<BC_dof_force> forces;

  // List of prescribed values for DoFs. This is old fashion realisation and will be changed soon.
  // List is populated by FEStorage::addBoundaryCondition(..). DoFs with fixed values are treated in nla3d in a
  // special manner: this DoFs are eliminated from global solve system and then reaction forces of a such DoFs
  // are found. For details see FEStorage::initializeSolutionData(..) and FEStorage::updateSolutionResults().
	std::list<BC_dof_constraint> constraints;

  // List of MPC equations. List is populated by FEStorage::addBoundaryCondition(..). An instance of Mpc class
  // is created outside of FEStorage class. But after add_bound(..) function FEStorage takes control on 
  // the instance of Mpc. MPCs can be deleted by FEStorage::deleteMpcs(). Of course, they are deleted in
  // ~FEStorage() destructor. 
	std::list<Mpc*> mpcs;

  // As in ansys mechanical APDL, here is a FE components - lists constist of number of entities. Here is the support
  // of nodal components and element components. This components can be used to apply BCs and/or MPCs. FE component 
  // is attached to FE Storage by FEStorage::addFEComponent(..) function. All stored components can be listed by
  // FEStorage::listFEComponents(..). One can access to the particular component by methods FEStorage:: getFEComponent(..). 
  std::vector<FEComponent*> feComponents;

  // An array of MpcCollection. Some number of MPCs can be grouped together because of they were produced 
  // by a single rule. For example, in FE model can rigid body MPC equations between the master node and slave nodes.
  // The relation between the master node and slave nodes described by a group of MPC equations. All of them are 
  // grouped in MpcCollection. See Mpc.h for more details. 
  // Here is a rule: Mpc defined inside of MpcCollection should be added to FEStorage::mpcs by add_boudns(..).
  // That means that Mpc instances are dynamically created in MpcCollection, but then FEStorage class takes control on it.
  // And FEStorage deletes Mpc by itself when it's needed.
  std::vector<MpcCollection*> mpcCollections;
  // TODO: here we store block matrix as 4 independent sparse matrices. It would be better to store
  // all of them in block matrix class (need to develop)
  // Matrix to be solved = [Kss,MPCsT;MPCs,0] - see details of how nla3d treats global system of equations 
  // in FEStorage::initializeSolutionData()
  math::SparseSymMatrix* KssMPCsT = nullptr; 
  // MPCc - a part of the global system of equations, MPC equation constants for constrained DoFs
  math::SparseMatrix* MPCc = nullptr; 
  // A part of the global K matrix with dimensions [numberOfConstrainedDofs x numberOfUnknownDofs]. It constaints of 
  // stiffness coefficients between constrained and unconstreined DoFs.
  math::SparseMatrix* Kcs = nullptr; 
  // A part of the global K matrix with dimensions [numberOfConstrainedDofs x numberOfConstrainedDofs]. It containts
  // of stiffness coefficients between constrained and constrained DoFs.
  math::SparseSymMatrix* Kcc = nullptr; 

  // matrices for transient analysis:
  // Damping matrix
  math::SparseSymMatrix* CssMPCsT = nullptr; 
  math::SparseMatrix* Ccs = nullptr; 
  math::SparseSymMatrix* Ccc = nullptr;
  // Damping matrix
  math::SparseSymMatrix* MssMPCsT = nullptr; 
  math::SparseMatrix* Mcs = nullptr; 
  math::SparseSymMatrix* Mcc = nullptr; 

  // A found values after converged loadstep. solutionValues has a size of numberOfDofs + numberOfMpcEq.
  // The whole vector is divided into different parts: solutionValues = [dofValues; mpcLagrangianValues],
  // where dofValues - FE model's DoFs, and mpcLagrangianValues - largrangian coefficients for solved MPC.
  // And dofValues is also represented by two parts: dofValues = [constrainedDofValues; unknownDofValues],
  // where constrainedDofValues - is a part for constrained DoFs, and unknownDofValues - is a part for solved dofs.
  // NOTE that the memory is allocated only once for solutionValues pointer. And others pointers 
  // just point to different starting addresses of solutionValues.
  std::vector<double> solutionValues;
  std::vector<double>::iterator dofValues;
  std::vector<double>::iterator constrainedDofValues;
  std::vector<double>::iterator unknownDofValues;
  std::vector<double>::iterator mpcLagrangianValues;

  // By incremetal approach a solution of a loadstep is found by several equilibrium iterations. 
  // solutionDeltaValues stores results (increment) of last equilibrium iteration. solutionDeltaValues,
  // dofDeltaValues, constrainedDofDeltaValues, unknownDofDeltaValues, mpcLagrangianDeltaValues
  // work the same as vectors described above (solutionValues,..)
  std::vector<double> solutionDeltaValues;
  std::vector<double>::iterator dofDeltaValues;
  std::vector<double>::iterator constrainedDofDeltaValues;
  std::vector<double>::iterator unknownDofDeltaValues;
  std::vector<double>::iterator mpcLagrangianDeltaValues;

  // Vector to store reactions for constrained DoFs. Size of [numberOfConstrainedDofs]. Reactions are calculated
  // after the global system solution in FEStorage::updateSolutionResults().
  std::vector<double> constrainedDofReactions;

  // vector externalForces stores external forces for Dofs. Size [numberOfDofs]. Pointers constrainedDofExternalForces,
  // unknownDofExternalForces just point to different parts of externalForces. External forces are defined
  // by boundary conditions in FEStorage::applyBoundaryConditions(..)
  std::vector<double> externalForces;
  std::vector<double>::iterator constrainedDofExternalForces;
  std::vector<double>::iterator unknownDofExternalForces;

  // Vector of constants for MPC equations. Size [numberOfMpcEq]
  std::vector<double> mpcConstantValues;

  // Vector rhs stores right hand side of the system of linear algebraic equations.
  // As long as constrained DoFs are eliminated from global system of equations, right
  // hand side of equations differ from unknownDofExternalForces and mpcConstantValues. Size of 
  // rhs [numberOfUnknownDofs+numberOfMpcEq]. See explanation in FEStorage::initializeSolutionData().
  // As above, unknownDofRhs, mpcEquationRhs are just pointers to different parts of rhs.
  std::vector<double> rhs;
  std::vector<double>::iterator unknownDofRhs;
  std::vector<double>::iterator mpcEquationRhs;

  // * In nla3d every DoF is represented by Dof class.
  // * DoF can be nodal or for element.
  // * elementDofs and nodeDofs are DofCollection 
  //    DofCollection keeps Dof objects and give it by (node, dof) pair
  // * Elements (or another entities) should register DoFs by using
  //    FEStorage::registerNodeDof(node, dof) or FEStorage::registerElementDof(element, dof).
  DofCollection elementDofs;
  DofCollection nodeDofs;

  // topology[node_number] = set of elements attached to it
  std::vector<std::set<uint32> > topology;

  // if transient is true that means that assembleGlobalEqMatrices should assemble M and C matrices
  // too
  bool transient;
};

// read Ansys Mechanical APDL *.cdb file. Nodes, Elements, Displacement BC and MPC (Constraint equations) is supported
// readCdbFile repcales storage's mesh.
bool readCdbFile (const char *filename, FEStorage *storage, ElementType elType);
} // namespace nla3d 

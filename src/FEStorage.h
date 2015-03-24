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
#include "PostProcessor.h"
#include "Mpc.h"
 
namespace nla3d {

//pre-defines
class Element;
class Node;
class PostProcessor;
class Dof;


#define ST_INIT 1
#define ST_LOADED 2
#define ST_SOLVED 3


// class FEStorage - heart of the NLA programm, it contains all data to keep in memory: 
//				nodes, elements, solution staff (global matrixes and vectors), processors.
//TODO: add mutexes to prevent multithreads errors

class FEStorage {
public:
	std::vector<Element*> elements;
	std::vector<Node> nodes;
	std::list<BC_dof_force> list_bc_dof_force;
	std::list<BC_dof_constraint> list_bc_dof_constraint;
	std::list<Mpc*> list_bc_MPC;
	std::vector<PostProcessor*> post_procs;
    std::vector<FEComponent*> feComponents;
    std::vector<MpcCollection*> mpcCollections;

	uint32 n_nodes;
	uint32 n_elements;
	uint32 n_dofs; //общее число степеней свободы системы

	uint32 n_constrained_dofs;
	uint32 n_solve_dofs;
	uint32 n_MPC_eqs;

	Material* material;
  nla3d::math::SparseSymmetricMatrix* KssCsT; //solve matrix [Kss,CsT;Cs,0]
  nla3d::math::SparseMatrix *Cc; // MPC eqs matrix for constrained dofs
  nla3d::math::SparseMatrix *Kcs; //part of global K matrix dim[n_constrained_dofs x n_solve_dofs]
  nla3d::math::SparseSymmetricMatrix *Kcc; //part of global K matrix dim[n_constrained_dofs x n_constrained_dofs]

	double* vec_q_lambda;
	double* vec_q;
	double* vec_qc;
	double* vec_qs;
	double* vec_lambda;

	double* vec_dq_dlambda;
	double* vec_dq;
	double* vec_dqc;
	double* vec_dqs;
	double* vec_dlambda;

	double* vec_reactions;

	double* vec_F;
	double* vec_Fc;
	double* vec_Fs;

	double* vec_b;

	double* vec_rhs;
	double* vec_rhs_dof;
	double* vec_rhs_mpc;

  // vector with size numberOfElements() * Element::numberOfDofs()
  std::vector<Dof*> elementDofs;
  // number of not NULL values in elementDofs
  uint32 numberOfElementDofs;
  // vector with size numberOfNodes() * Node::numberOfDofs()
  std::vector<Dof*> nodeDofs;
  // number of not NULL values in nodeDofs
  uint32 numberOfNodeDofs;

	uint16 status;

  ElementFactory::elTypes elType;

	//methods
	FEStorage();
	~FEStorage();

	uint32 getDofEqNumber(int32 node, Dof::dofType dof);
  Dof* getElementDof (int32 el, Dof::dofType dof);
  Dof* getNodeDof (int32 node, Dof::dofType dof);
  Dof* getDof (int32 node, Dof::dofType dof);
  bool isDofUsed (int32 node, Dof::dofType dof);
  void registerNodeDof(int32 node, Dof::dofType dof);
  void registerElementDof(int32 el, Dof::dofType dof);


  void Kij_add(int32 nodei, Dof::dofType dofi, int32 nodej, Dof::dofType dofj, double value);
  void Cij_add(uint32 eq_num, int32 nodej, Dof::dofType dofj, double coef);
  void Fi_add(int32 nodei, Dof::dofType dofi, double value);
	void zeroK();
	void zeroF();
  nla3d::math::SparseSymmetricMatrix& get_solve_mat();
	double* get_solve_rhs ();
	double* get_solve_result_vector();
	double* get_vec_dqs();
	Material& getMaterial();
	Node& getNode(uint32 _nn);
	Element& getElement(uint32 _en);
	PostProcessor& getPostProc(size_t _np);

	uint32 getNumberOfNodes () {
		return n_nodes;
	}
	uint32 getNumberOfElements () {
		return n_elements;
	}
	size_t getNumPostProc () {
		return post_procs.size();
	}
	uint16 getStatus() {
		return status;
	}
	void setStatus(uint16 st) {
		status = st;
	}
	uint32 get_n_solve_dofs() {
		return n_solve_dofs;
	}
	uint32 get_n_constrained_dofs() {
		return n_constrained_dofs;
	}
	uint32 gen_n_MPC_eqs() {
		return n_MPC_eqs;
	}

	std::list<BC_dof_constraint>& get_dof_const_BC_list();

	uint16 add_post_proc (PostProcessor *pp);
	void add_bounds (BC_dof_constraint &bc) ;
	void add_bounds (BC_dof_force &bc) ;
	void add_bounds (Mpc* bc);
  void addFEComponent (FEComponent* comp);
  void listFEComponents (std::ostream& out);
  FEComponent* getFEComponent(size_t i);
  FEComponent* getFEComponent(std::string name);

	void clearMesh ();
  void deleteElements();
	void nodes_reassign(uint32 nn); 
	void elements_reassign(uint32 en);

	bool prepare_for_solution ();

	void element_nodes(uint32 el, Node** node_ptr);
//	void get_q_e(uint32 el, double* ptr);
//	void get_dq_e(uint32 el, double* ptr);
  double getDofSolution (int32 node, Dof::dofType dof);
  double getDofSolutionDelta (int32 node, Dof::dofType dof);

	void getNodePosition(uint32 n, double* ptr, bool deformed = false); //def. or initial pos

	double getReaction(int32 n, Dof::dofType dof);

	void pre_first();
	void post_first();
	void process_solution();
	void apply_BCs (uint16 curLoadstep, uint16 curIteration, double d_par, double cum_par);
private:
};



inline double* FEStorage::get_solve_result_vector ()
{
	assert(vec_dq_dlambda);
	return vec_dqs;
}


inline Material& FEStorage::getMaterial()
{
	assert(material);
	return *material;
}

//getNode(nn)
//нумерация с 1
inline Node& FEStorage::getNode(uint32 _nn)
{
	assert(_nn <= n_nodes);
	return nodes[_nn-1];
}

//getElement(_en)
//нумерация с 1
inline Element& FEStorage::getElement(uint32 _en)
{
	assert(_en <= n_elements);
	return *(elements[_en-1]);
}

//getPostProc(_np)
//нумерация с 0
inline PostProcessor& FEStorage::getPostProc(size_t _np)
{
	assert(_np < getNumPostProc());
	return *post_procs[_np];
}


inline std::list<BC_dof_constraint>& FEStorage::get_dof_const_BC_list()
{
	return list_bc_dof_constraint;
}


//add_Bounds
inline void FEStorage::add_bounds (BC_dof_constraint &bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
	list_bc_dof_constraint.push_back(bc);
	n_constrained_dofs++;
}


inline void FEStorage::add_bounds (BC_dof_force &bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
	list_bc_dof_force.push_back(bc);
}


inline void FEStorage::add_bounds (Mpc* bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
    //TODO: we need to clear a memory after work 
	list_bc_MPC.push_back(bc);
}



inline double* FEStorage::get_vec_dqs()
{
	assert(vec_dq_dlambda);
	return vec_dqs;
}

bool read_ans_data (const char *filename, FEStorage *storage);

} // namespace nla3d 

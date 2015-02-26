#pragma once
#include <string>
#include <list>
#include "sys.h"
#include "Material.h"
#include "Element.h"
#include "math\Mat_Band.h"
#include "math\Sparse_Matrix.h"
#include "post_proc.h"
 
//pre-defines
class Element;
class Node;
class Post_proc;
class Dof;

class BC
{
public:
	BC () : is_updatetable(false)
	{  }
	//virtual void apply(FE_Storage_Interface *storage)=0;
	bool is_updatetable; //нужно ли обновлять на каждом шаге решения
};
class BC_dof_constraint : public BC
{
public:
	BC_dof_constraint () : node(0), node_dof(0), value(0.0)
	{	}
	int32 node;
	uint16 node_dof;
	double value;
	//void apply(FE_Storage_Interface *storage);
};

class BC_dof_force : public BC
{
public:
	BC_dof_force () : node(0), node_dof(0), value(0.0)
	{	}
	int32 node;
	uint16 node_dof;
	double value;
	//void apply(FE_Storage_Interface *storage);
};
class MPC_token
{
public:
	MPC_token() : node(0), node_dof(0), coef(0.0)
	{	}
	MPC_token(int32 n, uint16 dof, double coef) : node(n), node_dof(dof), coef(coef)
	{	}
	int32 node;
	uint16 node_dof;
	double coef;
};
// MPC equation like
// coef1 * dof1 + coef2 * dof2 + ... + coefn * dofn - b = 0
class BC_MPC : public BC
{
public:
	list<MPC_token> eq;
	double b;
	void apply(FE_Storage_Interface *storage);
};

#define ST_INIT 1
#define ST_LOADED 2
#define ST_SOLVED 3


// class FE_Storage - heart of the NLA programm, it contains all data to keep in memory: 
//				nodes, elements, solution staff (global matrixes and vectors), processors.
//TODO: add mutexes to prevent multithreads errors

class FE_Storage:
{
public:
	Element** elements;
	vector<Node> nodes;
	list<BC_dof_force> list_bc_dof_force;
	list<BC_dof_constraint> list_bc_dof_constraint;
	list<BC_MPC> list_bc_MPC;
	vector<Post_proc*> post_procs;

	uint32 n_nodes;
	uint32 n_elements;
	uint32 n_dofs; //общее число степеней свободы системы

	uint32 n_constrained_dofs;
	uint32 n_solve_dofs;
	uint32 n_MPC_eqs;

	Material *material;
	Sparse_Matrix_SUR *KssCsT; //solve matrix [Kss,CsT;Cs,0]
	Sparse_Matrix_R *Cc; // MPC eqs matrix for constrained dofs
	Sparse_Matrix_R *Kcs; //part of global K matrix dim[n_constrained_dofs x n_solve_dofs]
	Sparse_Matrix_SUR *Kcc; //part of global K matrix dim[n_constrained_dofs x n_constrained_dofs]

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

	Dof *dof_array; //массив степеней свободы, поизводит привязку номера степени свободы и номера уравнения
	//status
	uint16 status;

  Element::elTypes elType;

	//methods
	FE_Storage();
	~FE_Storage();

	uint32 get_dof_num(int32 node, uint16 dof);
	uint32 get_dof_eq_num(int32 node, uint16 dof);
	bool is_dof_constrained(int32 node, uint16 dof);

	void Kij_add(int32 nodei, uint16 dofi, int32 nodej, uint16 dofj, double value);
	void Cij_add(uint32 eq_num, int32 nodej, uint32 dofj, double value);
	void Fi_add(int32 nodei, uint16 dofi, double value);
	void zeroK();
	void zeroF();
	Sparse_Matrix_SUR& get_solve_mat();
	double* get_solve_rhs ();
	double* get_solve_result_vector();
	double* get_vec_dqs();
	Material& getMaterial();
	Node& getNode(uint32 _nn);
	Element& getElement(uint32 _en);
	Post_proc& getPostProc(uint32 _np);

	uint32 getNumDofs () {
		return n_dofs; 
	}
	uint32 getNumNode () {
		return n_nodes;
	}
	uint32 getNumElement () {
		return n_elements;
	}
	uint32 getNumPostProc () {
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

	list<BC_dof_constraint>& get_dof_const_BC_list();

	uint16 add_post_proc (Post_proc *pp);
	void add_bounds (BC_dof_constraint &bc) ;
	void add_bounds (BC_dof_force &bc) ;
	void add_bounds (BC_MPC &bc);

	void clearMesh ();
  void deleteElements();
	void nodes_reassign(uint32 nn); 
	void elements_reassign(uint32 en);

	bool prepare_for_solution ();

	void element_nodes(uint32 el, Node** node_ptr);
	void get_q_e(uint32 el, double* ptr);
	void get_dq_e(uint32 el, double* ptr);
	void get_q_n(uint32 n, double* ptr);
	double get_qi_n(int32 n, uint16 dof);
	void get_node_pos(uint32 n, double* ptr, bool def = false); //def. or initial pos

	double get_reaction_force(int32 n, uint16 dof);

	void pre_first();
	void post_first();
	void process_solution();
	void apply_BCs (uint16 curLoadstep, uint16 curIteration, double d_par, double cum_par);
private:
};


inline uint32 FE_Storage::get_dof_num(int32 node, uint16 dof) {
	// возвращает число от 1 до n_dofs
	uint32 res = (node < 0)?((-node-1)*Element::n_dofs()+dof+1):(n_elements*Element::n_dofs()+(node-1)*Node::n_dofs()+dof+1);
	return res; 
}

inline uint32 FE_Storage::get_dof_eq_num(int32 node, uint16 dof) {
	assert(dof_array);
	return dof_array[get_dof_num(node,dof)-1].eq_number;
}


inline bool FE_Storage::is_dof_constrained(int32 node, uint16 dof) {
	assert(dof_array);
	return dof_array[get_dof_num(node,dof)-1].is_constrained;
}


inline double* FE_Storage::get_solve_result_vector ()
{
	assert(vec_dq_dlambda);
	return vec_dqs;
}


inline Material& FE_Storage::getMaterial()
{
	assert(material);
	return *material;
}

//getNode(nn)
//нумерация с 1
inline Node& FE_Storage::getNode(uint32 _nn)
{
	assert(_nn <= n_nodes);
	return nodes[_nn-1];
}

//getElement(_en)
//нумерация с 1
inline Element& FE_Storage::getElement(uint32 _en)
{
	assert(_en <= n_elements);
	return *(elements[_en-1]);
}

//getPostProc(_np)
//нумерация с 0
inline Post_proc& FE_Storage::getPostProc(uint32 _np)
{
	assert(_np < getNumPostProc());
	return *post_procs[_np];
}


inline list<BC_dof_constraint>& FE_Storage::get_dof_const_BC_list()
{
	return list_bc_dof_constraint;
}


//add_Bounds
inline void FE_Storage::add_bounds (BC_dof_constraint &bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
	list_bc_dof_constraint.push_back(bc);
	n_constrained_dofs++;
}


inline void FE_Storage::add_bounds (BC_dof_force &bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
	list_bc_dof_force.push_back(bc);
}


inline void FE_Storage::add_bounds (BC_MPC &bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
	list_bc_MPC.push_back(bc);
}


// element_nodes(el, node_ptr), вызывающая сторона должна предоставить массив 
// указателей Node* на >= Element::nNodes() элементов
// el начинается с 1
inline void FE_Storage::element_nodes(uint32 el, Node** node_ptr)
{
	assert(el <= n_elements);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		node_ptr[i] = & nodes[elements[el-1]->node_num(i)-1]; //TODO: CHECK
}


//get_qi_n(n, dof)
// n начинается с 1
inline double FE_Storage::get_qi_n(int32 n, uint16 dof)
{
	assert(vec_q_lambda);
	return vec_q[get_dof_eq_num(n, dof)-1];
}


inline double* FE_Storage::get_vec_dqs()
{
	assert(vec_dq_dlambda);
	return vec_dqs;
}

bool read_ans_data (const char *filename, FE_Storage_Interface *storage);

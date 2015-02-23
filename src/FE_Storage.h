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
//struct Bound 
//{
//	uint32 node;
//	uint16 key;
//	double value;
//	void display(uint32 bn)
//	{
//		string k="UNKNOWN";
//		switch (key)
//		{
//			case D_UX: k = "UX";
//				break;
//			case D_UY: k = "UY";
//				break;
//			case D_UZ: k = "UZ";
//				break;
//		}
//		echo("B %d: %d\t%s\t%f",bn, node, k.c_str(),value);
//	}
//	string toString()
//	{
//		string str;
//		char buff[100];
//		switch (key)
//		{
//			case D_UX: str="UX";
//				break;
//			case D_UY: str="UY";
//				 break;
//			case D_UZ: str="UZ";
//				 break;
//			default:
//				warning("Bound::toString: unknown bound key %d", key);
//		}
//		sprintf_s(buff,100," %d %f",node, value);
//		str+=buff;
//		return str;
//	}
//	void read_from_stream (istream &file)
//	{
//		string str;
//		file >> str;
//		if (str=="UX")
//			key = D_UX;
//		else if (str=="UY")
//			key = D_UY;
//		else if (str=="UZ")
//			key = D_UZ;
//		else
//		{
//			//TODO: error
//			warning("Bound::read_from_stream: unknown bound key %d", key);
//		}
//		file >> node;
//		file >> value;
//	}
//};

#define ST_INIT 1
#define ST_LOADED 2
#define ST_SOLVED 3

// abstract class to provide interface to Storage
class FE_Storage_Interface
{
public:
	virtual uint32 get_dof_num(int32 node, uint16 dof) = 0;
	virtual uint32 get_dof_eq_num(int32 node, uint16 dof) = 0;
	virtual bool is_dof_constrained(int32 node, uint16 dof) = 0;
	virtual void Kij_add(int32 nodei, uint16 dofi, int32 nodej, uint16 dofj, double value) = 0;
	virtual void Cij_add(uint32 eq_num, int32 nodej, uint32 dofj, double value) = 0;
	virtual void Fi_add(int32 nodei, uint16 dofi, double value) = 0;
	virtual void zeroK() = 0;
	virtual void zeroF() = 0;
	virtual	Sparse_Matrix_SUR& get_solve_mat() = 0;
	virtual double* get_solve_rhs () = 0;
	virtual double* get_solve_result_vector() = 0;
	virtual double* get_vec_dqs() = 0;

	virtual Material& getMaterial() = 0;
	virtual Node& getNode(uint32 _nn) = 0;
	virtual Element& getElement(uint32 _en) = 0;
	virtual Post_proc& getPostProc(uint32 _np) = 0;

	virtual uint32 get_nBand () = 0;
	virtual uint32 getNumDofs () = 0;
	virtual uint32 getNumNode () = 0;
	virtual uint32 getNumElement () = 0;
	virtual uint32 getNumPostProc () = 0;
	virtual uint32 get_n_solve_dofs() = 0;
	virtual uint32 get_n_constrained_dofs() = 0;
	virtual uint32 gen_n_MPC_eqs() = 0;
	virtual uint16 getStatus() = 0;
	//virtual double get_dof_solution(int32 node, uint16 dof);
	virtual list<BC_dof_constraint>& get_dof_const_BC_list()=0;
	virtual void setStatus(uint16 st) = 0;

	virtual uint16 add_post_proc (Post_proc *pp) = 0;
	virtual void add_bounds (BC_dof_constraint &bc) = 0;
	virtual void add_bounds (BC_dof_force &bc) = 0;
	virtual void add_bounds (BC_MPC &bc) = 0;

	virtual void clearMesh () = 0;
	virtual void nodes_reassign(uint32 nn) = 0; 
	virtual void elements_reassign(uint32 en) = 0;

	virtual bool prepare_for_solution ()  =0;

	virtual void element_nodes(uint32 el, Node** node_ptr) = 0;
	virtual void get_q_e(uint32 el, double* ptr) = 0;
	virtual void get_dq_e(uint32 el, double* ptr) = 0;
	virtual void get_q_n(uint32 n, double* ptr) = 0;
	virtual double get_qi_n(int32 n, uint16 dof) = 0;
	virtual void get_node_pos(uint32 n, double* ptr, bool def = false) = 0; //def. or initial pos
	virtual double get_reaction_force(int32 n, uint16 dof)=0;


	virtual void pre_first() = 0;
	virtual void post_first() = 0;
	virtual void process_solution() = 0;
	virtual void apply_BCs (uint16 curLoadstep, uint16 curIteration, double d_par, double cum_par) = 0;
};

// class FE_Storage - heart of the NLA programm, it contains all data to keep in memory: 
//				nodes, elements, solution staff (global matrixes and vectors), processors.
//TODO: add mutexes to prevent multithreads errors
template <typename el_type>
class FE_Storage : public FE_Storage_Interface
{
public:
	FE_Storage();
	~FE_Storage();
	el_type dummy; // пустышка, для разных нужд static initialization
	vector<el_type> elements;
	vector<Node> nodes;
	//vector<Bound> bounds;
	list<BC_dof_force> list_bc_dof_force;
	list<BC_dof_constraint> list_bc_dof_constraint;
	list<BC_MPC> list_bc_MPC;
	vector<Post_proc*> post_procs;

	uint32 n_nodes;
	uint32 n_elements;
	uint32 n_dofs; //общее число степеней свободы системы
	uint32 nBand;

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

	//реализация интерфейса FE_Storage_Interface

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

	uint32 get_nBand ();
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
	void count_nBand ();
};


template <typename el_type>
FE_Storage<el_type>::FE_Storage() : dummy()
{
	n_nodes = 0;
	n_elements = 0;
	n_dofs = 0;
	nBand = 0;

	n_constrained_dofs = 0;
	n_solve_dofs = 0;
	n_MPC_eqs = 0;

	material = NULL;
	KssCsT = NULL;
	Cc = NULL;
	Kcs = NULL;
	Kcc = NULL;

	vec_q_lambda = NULL;
	vec_q = NULL;
	vec_qc = NULL;
	vec_qs = NULL;
	vec_lambda = NULL;

	vec_dq_dlambda = NULL;
	vec_dq = NULL;
	vec_dqc = NULL;
	vec_dqs = NULL;
	vec_dlambda = NULL;

	vec_reactions = NULL;

	vec_F = NULL;
	vec_Fc = NULL;
	vec_Fs = NULL;

	vec_b = NULL;

	vec_rhs = NULL;
	vec_rhs_dof = NULL;
	vec_rhs_mpc = NULL;

	dof_array = NULL;
	
	status = ST_INIT;
};
template <typename el_type>
FE_Storage<el_type>::~FE_Storage ()
{
	status = ST_INIT; //TODO: check
	for (uint16 i=0; i < post_procs.size(); i++)
		delete post_procs[i];
	if (material) delete material;
	if (KssCsT) delete KssCsT;
	if (Cc) delete Cc;
	if (Kcs) delete Kcs;
	if (Kcc) delete Kcc;

	if (vec_q_lambda) delete[] vec_q_lambda;
	if (vec_dq_dlambda) delete[] vec_dq_dlambda;
	if (vec_reactions) delete[] vec_reactions;
	if (vec_F) delete[] vec_F;
	if (vec_b) delete[] vec_b;
	if (vec_rhs) delete[] vec_rhs;
	if (dof_array) delete[] dof_array;

	elements.clear();
	nodes.clear();
	list_bc_dof_force.clear();
	list_bc_dof_constraint.clear();
	list_bc_MPC.clear();
}

template <typename el_type>
inline uint32 FE_Storage<el_type>::get_dof_num(int32 node, uint16 dof)
{
	// возвращает число от 1 до n_dofs
	uint32 res = (node < 0)?((-node-1)*Element::n_dofs()+dof+1):(n_elements*Element::n_dofs()+(node-1)*Node::n_dofs()+dof+1);
	return res; 
}

template <typename el_type>
inline uint32 FE_Storage<el_type>::get_dof_eq_num(int32 node, uint16 dof)
{
	assert(dof_array);
	return dof_array[get_dof_num(node,dof)-1].eq_number;
}

template <typename el_type>
inline bool FE_Storage<el_type>::is_dof_constrained(int32 node, uint16 dof)
{
	assert(dof_array);
	return dof_array[get_dof_num(node,dof)-1].is_constrained;
}

template <typename el_type>
bool FE_Storage<el_type>::prepare_for_solution ()
{
	if (status < ST_LOADED)
	{
		warning ("FE_Storage::prepare_for_solution: model isn't loaded");
		return false;
	}

	n_dofs = n_nodes*Node::n_dofs() + n_elements*Element::n_dofs(); //общее число степеней свободы
	//n_constrained_dofs - подсчитано при задании BC's
	//n_MPC_eqs - подсчитано при задании MPC's
	n_MPC_eqs = list_bc_MPC.size();
	//
	n_solve_dofs = n_dofs - n_constrained_dofs;
	if (n_MPC_eqs)	Cc = new Sparse_Matrix_R(n_MPC_eqs, n_constrained_dofs);
	KssCsT = new Sparse_Matrix_SUR(n_solve_dofs + n_MPC_eqs);
	Kcs = new Sparse_Matrix_R(n_constrained_dofs, n_solve_dofs);
	Kcc = new Sparse_Matrix_SUR(n_constrained_dofs);

	// заполняем массив степеней свободы
	dof_array = new Dof[n_dofs];
	uint32 next_eq_solve = n_constrained_dofs+1;
	uint32 next_eq_const = 1;
	list<BC_dof_constraint>::iterator p = list_bc_dof_constraint.begin();
	while (p != list_bc_dof_constraint.end())
	{
		dof_array[get_dof_num(p->node, p->node_dof)-1].is_constrained = true;
		p++;
	}
	for (uint32 i=0; i<n_dofs; i++)
	{
		if (dof_array[i].is_constrained)
			dof_array[i].eq_number = next_eq_const++;
		else
			dof_array[i].eq_number = next_eq_solve++;
	}
	//DEBUG: print dof array
	//cout << "Dof_array:" << endl;
	//for (uint32 i=0; i<n_dofs; i++)
	//{
	//	cout << i+1 << ": " << dof_array[i].eq_number << ", " << dof_array[i].is_constrained << endl;
	//}
	//настраеваем массивы значенийстепеней свободы {qc; qs; lambda}
	vec_q_lambda = new double[n_dofs+n_MPC_eqs];
	memset(vec_q_lambda, 0, sizeof(double)*(n_dofs+n_MPC_eqs));
	vec_q = vec_q_lambda;
	vec_qc = vec_q_lambda;
	vec_qs = &(vec_q_lambda[n_constrained_dofs]);
	vec_lambda = &(vec_q_lambda[n_dofs]); //тут нарушается адресация
	
	vec_dq_dlambda = new double[n_dofs+n_MPC_eqs];
	memset(vec_dq_dlambda, 0, sizeof(double)*(n_dofs+n_MPC_eqs));
	vec_dq = vec_dq_dlambda;
	vec_dqc = vec_dq_dlambda;
	vec_dqs = &(vec_dq_dlambda[n_constrained_dofs]);
	vec_dlambda = &(vec_dq_dlambda[n_dofs]);
	
	vec_reactions = new double[n_constrained_dofs];
	memset(vec_reactions, 0, sizeof(double)*n_constrained_dofs);

	vec_F = new double[n_dofs];
	memset(vec_F, 0, sizeof(double)*n_dofs);
	vec_Fc = vec_F;
	vec_Fs = &(vec_F[n_constrained_dofs]);

	if (n_MPC_eqs)
	{
		vec_b = new double[n_MPC_eqs];
		memset(vec_b, 0, sizeof(double)*n_MPC_eqs);
	}

	// {rhs} = size(n_solve_dofs+n_MPC_eqs)
	vec_rhs = new double[n_solve_dofs+n_MPC_eqs];
	memset(vec_rhs, 0, sizeof(double)*(n_solve_dofs+n_MPC_eqs));
	vec_rhs_dof = vec_rhs;
	vec_rhs_mpc = &(vec_rhs[n_solve_dofs]);

	echolog("DoFs = %d, constrained DoFs = %d, MPC eq. = %d, TOTAL eq. = %d", n_dofs, n_constrained_dofs, n_MPC_eqs, n_solve_dofs + n_MPC_eqs);

	if (!material)
	{
		warning("FE_Storage::prepare_for_solution: material isn't defined");
		return false;
	}
	return true;
}

template <typename el_type>
uint32 FE_Storage<el_type>::get_nBand ()
{
	if (!nBand) count_nBand();
	return nBand;
}
//
template <typename el_type>
void FE_Storage<el_type>::count_nBand ()
{
	//TODO: после введения степеней свободы элемента эта функция работает непарвильно
	uint32 min, max, MAX;
	MAX = 0;
	for (uint32 i=0; i<n_elements; i++)
	{
		
		min = elements[i].node_num(0);
		max = min;
		for (uint16 j=1; j < Element::n_nodes(); j++) 
		{
			if (elements[i].node_num(j) < min)
				min = elements[i].node_num(j);
			if (elements[i].node_num(j) > max)
				max = elements[i].node_num(j);
		}
		if ((max-min) > MAX)
			MAX = max - min;
	}
	nBand = (MAX+1)*Node::n_dofs(); // полуширина ленты с центральной диагональю
}

//функция возвращает ссылку на значение, которое хранится в глобальной матрице жосткости по адресу 
// строка: dofi-тая степень свободы nodei-того узла 
// столбец: dofj-тая степень свободы nodej-того узла 
// узлы нумеруются с индекса 1, dofi=[0;Node::nDOF()-1];
template<typename el_type>
void FE_Storage<el_type>::Kij_add(int32 nodei, uint16 dofi, int32 nodej, uint16 dofj, double value)
{
	//assert(matK);
	uint32 eq_row = get_dof_eq_num(nodei, dofi);
	uint32 eq_col = get_dof_eq_num(nodej, dofj);
	if (eq_row > eq_col) swap(eq_row, eq_col);
	if (eq_row <= n_constrained_dofs)
		if (eq_col <= n_constrained_dofs)
			Kcc->add_value(eq_row, eq_col, value);
		else
			Kcs->add_value(eq_row, eq_col - n_constrained_dofs, value);
	else
		KssCsT->add_value(eq_row - n_constrained_dofs, eq_col - n_constrained_dofs, value);

}
template<typename el_type>
void FE_Storage<el_type>::Cij_add(uint32 eq_num, int32 nodej, uint32 dofj, double value)
{
	assert(Cc);
	assert(eq_num > 0 && eq_num <= n_MPC_eqs);
	uint32 dof_col = get_dof_eq_num(nodej, dofj);
	if (dof_col <= n_constrained_dofs) 
		Cc->add_value(eq_num, dof_col, value);
	else
		//Cs
		KssCsT->add_value(dof_col-n_constrained_dofs, n_solve_dofs+eq_num, value);
}
//getFi(..) возвращает ссылку на ячейку глобального вектора сил
//по адресу dofi-тая степень свободы nodei-того узла
// узлы нумеруются с индекса 1, dofi=[0;Node::nDOF()-1];
// TODO: CHECK
template<typename el_type>
void FE_Storage<el_type>::Fi_add(int32 nodei, uint16 dofi, double value)
{
	assert(vec_F);
	uint32 row = get_dof_eq_num(nodei, dofi);
	assert(row <= n_dofs);
	vec_F[row-1] += value;
}
// clearK() очищает матрицу жесткости
template<typename el_type>
void FE_Storage<el_type>::zeroK()
{
	assert(KssCsT);
	KssCsT->zero_block(n_solve_dofs);// TODO: сделать эту функцию //DEBUG
	//KssCsT->zero();//DEBUG
	Kcc->zero();
	Kcs->zero();
}
// clearF() очищает глобальный вектор сил
template<typename el_type>
void FE_Storage<el_type>::zeroF()
{
	assert(vec_F);
	memset(vec_F, 0, sizeof(double)*n_dofs);
	memset(vec_dq_dlambda, 0, sizeof(double)*(n_dofs+n_MPC_eqs));//TODO: test
}
// getMatK возвращает ссылку на структуру данных глобальной матрицы жесткости
template<typename el_type>
Sparse_Matrix_SUR& FE_Storage<el_type>::get_solve_mat()
{
	assert(KssCsT);
	return *KssCsT;
}

//getVecF()
template<typename el_type>
double* FE_Storage<el_type>::get_solve_rhs ()
{
	assert(vec_rhs);
	//собрать его тут
	//cout << endl << "vec_dqc" << endl; //DEBUG
	//for (int i=0; i < n_constrained_dofs; i++)
	//	cout << vec_dqc[i] << endl;
	
	//cout << endl << "vec_F" << endl; //DEBUG
	//for (int i=0; i < n_dofs; i++)
	//	cout << vec_F[i] << endl;
	double *KcsTdqc = new double[n_solve_dofs];
	Kcs->transpose_mult_vec(vec_dqc,KcsTdqc);
	for (uint32 i=0; i < n_solve_dofs; i++)
		vec_rhs_dof[i] = vec_Fs[i] - KcsTdqc[i];//Kcs->transpose_mult_vec_i(vec_dqc,i+1);//TODO: тут может быть ошибка
	for (uint32 i=0; i < n_MPC_eqs; i++)
		vec_rhs_mpc[i] = vec_b[i] - Cc->mult_vec_i(vec_dqc,i+1);
	delete[] KcsTdqc;
	//cout << endl << "vec_rhs_dof" << endl; //DEBUG
	//for (int i=0; i < n_solve_dofs; i++)
	//	cout << vec_rhs[i] << endl;

	//Kcc->print(); //DEBUG
	//KssCsT->print();
	//Kcs->print();
	return vec_rhs;
}


template<typename el_type>
double* FE_Storage<el_type>::get_solve_result_vector ()
{
	assert(vec_dq_dlambda);
	return vec_dqs;
}

//getMaterial()
template<typename el_type>
Material& FE_Storage<el_type>::getMaterial()
{
	assert(material);
	return *material;
}

//getNode(nn)
//нумерация с 1
template<typename el_type>
Node& FE_Storage<el_type>::getNode(uint32 _nn)
{
	assert(_nn <= n_nodes);
	return nodes[_nn-1];
}

//getElement(_en)
//нумерация с 1
template<typename el_type>
Element& FE_Storage<el_type>::getElement(uint32 _en)
{
	assert(_en <= n_elements);
	return elements[_en-1];
}

//getPostProc(_np)
//нумерация с 0
template<typename el_type>
Post_proc& FE_Storage<el_type>::getPostProc(uint32 _np)
{
	assert(_np < getNumPostProc());
	return *post_procs[_np];
}

//getBound(_nb)
//нумерация с 0
//template<typename el_type>
//Bound& FE_Storage<el_type>::getBound (uint32 _nb)
//{
//	assert(_nb < getNumBound());
//	return bounds[_nb];
//}

//add_post_proc(*pp)

template<typename el_type>
list<BC_dof_constraint>& FE_Storage<el_type>::get_dof_const_BC_list()
{
	return list_bc_dof_constraint;
}

template<typename el_type>
uint16 FE_Storage<el_type>::add_post_proc (Post_proc *pp)
{
	assert(pp);
	uint16 num = this->post_procs.size()+1;
	pp->nPost_proc = num;
	post_procs.push_back(pp);
	return num;
}

//add_Bounds
template<typename el_type>
void FE_Storage<el_type>::add_bounds (BC_dof_constraint &bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
	list_bc_dof_constraint.push_back(bc);
	n_constrained_dofs++;
}

template<typename el_type>
void FE_Storage<el_type>::add_bounds (BC_dof_force &bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
	list_bc_dof_force.push_back(bc);
}

template<typename el_type>
void FE_Storage<el_type>::add_bounds (BC_MPC &bc)
{
	//TODO: предполагается, что закрепления не повторяют одну и туже степень свободы!
	list_bc_MPC.push_back(bc);
}

//clearMesh() функция очищает таблицу узлов, элементов, ГУ
template<typename el_type>
void FE_Storage<el_type>::clearMesh ()
{
	status = ST_INIT; //TODO продумать
	elements.clear();
	nodes.clear();
	list_bc_dof_constraint.clear();
	list_bc_dof_force.clear();
	list_bc_MPC.clear();
	n_nodes = 0;
	n_elements = 0;
	n_dofs = 0;
	n_constrained_dofs = 0;
	n_solve_dofs = 0;
	n_MPC_eqs = 0;
	nBand = 0;
}

//nodes_reassign(_nn)
template<typename el_type>
void FE_Storage<el_type>::nodes_reassign(uint32 _nn)
{
	nodes.clear();
	n_nodes = _nn;
  //Node() fires Vec<3> constructor, thus Node coordinates are (0,0,0) by default
  //TODO: try-catch of memory overflow
	nodes.assign(_nn, Node()); // TODO: тут бы catch на возможность выделения памяти
}

//elements_reassign(_en)
template<typename el_type>
void FE_Storage<el_type>::elements_reassign(uint32 _en)
{
	elements.clear();
	n_elements = _en;
	elements.assign(_en, el_type());
	nBand = 0;
}

// element_nodes(el, node_ptr), вызывающая сторона должна предоставить массив 
// указателей Node* на >= Element::nNodes() элементов
// el начинается с 1
template<typename el_type>
void FE_Storage<el_type>::element_nodes(uint32 el, Node** node_ptr)
{
	assert(el <= n_elements);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		node_ptr[i] = & nodes[elements[el-1].node_num(i)-1]; //TODO: CHECK
}

// get_q_e(el, ptr) функция возвращает вектор узловых степеней свободы элемента,
// вызывающая сторона должна предоставить массив ptr размерностью Element::n_nodes()*Node::n_dofs() + Element::n_dofs()
// el начинается с 1
template<typename el_type>
void FE_Storage<el_type>::get_q_e(uint32 el, double* ptr)
{
	assert(el <= n_elements);
	assert(vec_q_lambda);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		for (uint16 j=0; j<Node::n_dofs(); j++)
			ptr[i*Node::n_dofs()+j] = vec_q_lambda[get_dof_eq_num(elements[el-1].node_num(i), j)-1];
	for (uint16 i=0; i < Element::n_dofs(); i++)
		ptr[Element::n_nodes()*Node::n_dofs()+i] = vec_q_lambda[get_dof_eq_num(-(int32)el, i)-1];
}
// get_dq_e см. выше
template<typename el_type>
void FE_Storage<el_type>::get_dq_e(uint32 el, double* ptr)
{
	assert(el <= n_elements);
	assert(vec_dq_dlambda);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		for (uint16 j=0; j<Node::n_dofs(); j++)
			ptr[i*Node::n_dofs()+j] = vec_dq_dlambda[get_dof_eq_num(elements[el-1].node_num(i), j)-1];
	for (uint16 i=0; i < Element::n_dofs(); i++)
		ptr[Element::n_nodes()*Node::n_dofs()+i] = vec_dq_dlambda[get_dof_eq_num(-(int32)el, i)-1];
}

//get_q_n(n, ptr)
// n начинается с 1
template<typename el_type>
void FE_Storage<el_type>::get_q_n(uint32 n, double* ptr)
{
	assert(n > 0 && n <= n_nodes);
	assert(vec_q_lambda);
	for (uint16 j=0; j<Node::n_dofs(); j++)
		ptr[j] = vec_q[get_dof_eq_num(n, j)-1];
}

//get_qi_n(n, dof)
// n начинается с 1
template<typename el_type>
double FE_Storage<el_type>::get_qi_n(int32 n, uint16 dof)
{
	assert(vec_q_lambda);
	return vec_q[get_dof_eq_num(n, dof)-1];
}

//void get_node_pos(uint32 n, double* ptr, bool def = false) double массим на 3 элемента!
//n с 1
template<typename el_type>
void FE_Storage<el_type>::get_node_pos(uint32 n, double* ptr, bool def = false)
{
	assert(n > 0 && n <= n_nodes);
	for (uint16 i=0; i<3; i++)
		ptr[i] = nodes[n-1].pos[i];
	if (def)
	{
		for (uint16 i=0; i<Node::n_dofs(); i++)
			switch(Node::dof_type(i))
			{
			case UX:
				ptr[0] += get_qi_n(n, i);
				break;
			case UY:
				ptr[1] += get_qi_n(n,i);
				break;
			case UZ:
				ptr[2] += get_qi_n(n,i);
			}
	}
}

template <typename el_type>
double FE_Storage<el_type>::get_reaction_force(int32 n, uint16 dof)
{
	uint32 eq_num = get_dof_eq_num(n,dof);
	assert(vec_reactions);
	assert(eq_num > 0 && eq_num <= n_constrained_dofs);
	return vec_reactions[eq_num-1];
}

template <typename el_type>
void FE_Storage<el_type>::pre_first()
{
	//включаем тренировку матриц
	KssCsT->start_training();
	if (n_MPC_eqs) Cc->start_training();
	Kcc->start_training();
	Kcs->start_training();
}

template <typename el_type>
void FE_Storage<el_type>::post_first()
{
	//выключаем тренировку матриц
	if (n_MPC_eqs) 
	{
		Cc->stop_training();
	}
	KssCsT->stop_training();
	
	Kcc->stop_training();
	Kcs->stop_training();

}

template <typename el_type>
void FE_Storage<el_type>::process_solution()
{
	// из вектора решения вытаскиваем все необходимое
	//складываем решение

	for (uint32 i=0; i < n_dofs; i++)
		vec_q[i] += vec_dq[i];
	for (uint32 i=n_dofs; i < n_dofs+n_MPC_eqs; i++)
		vec_q_lambda[i] = vec_dq_dlambda[i];
	
	//cout << "q:"<<endl;
	//for (uint32 i=0; i < n_nodes; i++)
	//	for (uint16 j=0; j < Node::n_dofs(); j++)
	//		cout << "N" <<i+1<<"("<<j<<")=" << get_qi_n(i+1, j) << endl; //DEBUG
	//находим реакции
	//cout << "reactions:"<<endl;
	for (uint32 i=0; i < n_constrained_dofs; i++)
	{
		vec_reactions[i] = Kcs->mult_vec_i(vec_dqs,i+1) + Kcc->mult_vec_i(vec_dqc,i+1) - vec_Fc[i];
		if (n_MPC_eqs && Cc->get_n_values())
			vec_reactions[i] += Cc->transpose_mult_vec_i(vec_dlambda,i+1);
		//cout << vec_reactions[i] << endl; //DEBUG
	}
}

template <typename el_type>
void FE_Storage<el_type>::apply_BCs (uint16 curLoadstep, uint16 curIteration, double d_par, double cum_par)
{
	if (curLoadstep == 1 && curIteration == 1)
	{
		//заполним Cc Cs и vec_b
		uint32 eq_num = 1;
		list<BC_MPC>::iterator mpc = list_bc_MPC.begin();
		while (mpc != list_bc_MPC.end())
		{
			vec_b[eq_num-1] = mpc->b;
			list<MPC_token>::iterator token = mpc->eq.begin();
			while (token != mpc->eq.end())
			{
				Cij_add(eq_num, token->node, token->node_dof, token->coef);
				token++;
			}
			eq_num++;
			mpc++;
		}
	}
	
	//заполним вектор узловых сил
	list<BC_dof_force>::iterator bc_force = list_bc_dof_force.begin();
	while (bc_force != list_bc_dof_force.end())
	{
		Fi_add(bc_force->node, bc_force->node_dof, bc_force->value*cum_par);
		bc_force++;
	}

	//заполним вектор заданных перемещений
	list<BC_dof_constraint>::iterator bc_dof = list_bc_dof_constraint.begin();
	while (bc_dof != list_bc_dof_constraint.end())
	{
		uint32 eq_num = get_dof_eq_num(bc_dof->node, bc_dof->node_dof);//TODO: DEBUG stuff
		vec_dq[eq_num-1] = bc_dof->value*d_par;
		bc_dof++;
	}
}
template <typename el_type>
double* FE_Storage<el_type>::get_vec_dqs()
{
	assert(vec_dq_dlambda);
	return vec_dqs;
}

bool read_ans_data (const char *filename, FE_Storage_Interface *storage);

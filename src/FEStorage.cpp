// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "FEStorage.h"

namespace nla3d {
using namespace math;

FEStorage::FEStorage()  {
	n_nodes = 0;
	n_elements = 0;
	n_dofs = 0;

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

	status = ST_INIT;
  elType = ElementFactory::NOT_DEFINED;

  numberOfElementDofs = 0;
  numberOfNodeDofs = 0;
};

FEStorage::~FEStorage () {
	status = ST_INIT;
	for (size_t i = 0; i < post_procs.size(); i++) {
		delete post_procs[i];
  }
  post_procs.clear();
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

  for (size_t i = 0; i < feComponents.size(); i++) {
    delete feComponents[i];
  }
  feComponents.clear();
  clearMesh();

  for (size_t i = 0; i < elementDofs.size(); i++) {
    if (elementDofs[i] != NULL) {
      delete elementDofs[i];
    }
  }
  elementDofs.clear();

  for (size_t i = 0; i < nodeDofs.size(); i++) {
    if (nodeDofs[i] != NULL) {
      delete nodeDofs[i];
    }
  }
  nodeDofs.clear();

  //TODO: we need to clear a memory after work 
  // list_bc_MPC
  // mpcCollections
}

bool FEStorage::prepare_for_solution () {
  assert(n_nodes > 0 && n_elements > 0);
	if (status < ST_LOADED) {
		warning ("FEStorage::prepare_for_solution: model isn't loaded");
		return false;
	}
  
  for (size_t i = 0; i < mpcCollections.size(); i++) {
    mpcCollections[i]->pre();
    mpcCollections[i]->registerMpcsInStorage();
  }

  // Total number of dofs (only registered by elements)
	n_dofs = numberOfElementDofs + numberOfNodeDofs;

  //TODO: make it possible to work without any constrained dofs
  // (in case of only MPC)
  n_constrained_dofs = static_cast<uint32> (list_bc_dof_constraint.size());
	n_MPC_eqs = static_cast<uint32> (list_bc_MPC.size());

  // n_solve_dofs - number of Dof need to be found on every step
	n_solve_dofs = n_dofs - n_constrained_dofs;
  // In nla3d solution procedure there are 3 distinguish types of unknowns. First one "c" - constrained degress of freedom,
  // and consequently known at solution time. Second one "s" - degrees of freedom need to be found (solved).
  // And last one "lambda" - lagrange multipliers for applied MPC constraints.
  // Cc [n_MPC_eqs x n_constrained_dofs]  - part of a global stiffness matrix with MPC coefficients for constrained dofs
  // A global System of Linera Algebraic Equations:
  //
  //  |  Kcc  |  Kcs  | Cc^T |   |   qc  |   |vec_Fc|
  //  |-------|--------------|   |-------|   |      |
  //  | Kcs^T |              | * |   qs  | = |vec_Fs|
  //  |-------|    KssCsT    |   |-------|   |      |
  //  |  Cc   |              |   | lambda|   |vec_b |
  //  
  //  1. But really only
  //
  //  KssCsT * [qs; lambda]^T = [vec_rhs_dof; vec_rhs_mpc]^T
  //
  //  to be solved by eq. solver. 
  //
  //  where:
  //  vec_rhs_dof = vec_Fs - Kcs^T*qc
  //  vec_rhs_mpc = vec_b - Cc*qc
  //
  //  2. Then to restore reaction forces:
  //  vec_Fc = Kcc*qc + Kcs*qs + Cc^T*lambda
  //
	if (n_MPC_eqs) {
    Cc = new math::SparseMatrix(n_MPC_eqs, n_constrained_dofs);
  }
	KssCsT = new math::SparseSymmetricMatrix(n_solve_dofs + n_MPC_eqs);
	Kcs = new math::SparseMatrix(n_constrained_dofs, n_solve_dofs);
	Kcc = new math::SparseSymmetricMatrix(n_constrained_dofs);

	// Fill elementsDofs and nodeDofs with isConstrained information
  // and then give them an equation numbers
	uint32 next_eq_solve = n_constrained_dofs+1;
	uint32 next_eq_const = 1;
	list<BC_dof_constraint>::iterator p = list_bc_dof_constraint.begin();
	while (p != list_bc_dof_constraint.end()) {
    getDof(p->node, p->node_dof)->isConstrained = true;
		p++;
	}
	for (size_t i = 0; i < getNumberOfElements()*Element::getNumberOfDofs(); i++)	{
    if (elementDofs[i] == NULL) {
      continue;
    }
		if (elementDofs[i]->isConstrained) {
			elementDofs[i]->eqNumber = next_eq_const++;
    } else {
			elementDofs[i]->eqNumber = next_eq_solve++;
    }
	}

	for (size_t i = 0; i < getNumberOfNodes()*Node::getNumberOfDofs(); i++)	{
    if (nodeDofs[i] == NULL) {
      continue;
    }
		if (nodeDofs[i]->isConstrained) {
			nodeDofs[i]->eqNumber = next_eq_const++;
    } else {
			nodeDofs[i]->eqNumber = next_eq_solve++;
    }
	}
  assert(next_eq_const - 1 == n_constrained_dofs);
  assert(next_eq_solve - 1 == n_dofs);

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

	if (n_MPC_eqs) {
		vec_b = new double[n_MPC_eqs];
		memset(vec_b, 0, sizeof(double)*n_MPC_eqs);
	}

	// {rhs} = size(n_solve_dofs+n_MPC_eqs)
	vec_rhs = new double[n_solve_dofs+n_MPC_eqs];
	memset(vec_rhs, 0, sizeof(double)*(n_solve_dofs+n_MPC_eqs));
	vec_rhs_dof = vec_rhs;
	vec_rhs_mpc = &(vec_rhs[n_solve_dofs]);

  std::stringstream ss;
  ss << "Nodal DoF types:";
  for (uint16 i = 0; i < Node::getNumberOfDofs(); i++) {
    ss << " " << Dof::dofTypeLabels[Node::getDofType(i)];
  }
  ss << " Number of nodal DoFs: " << numberOfNodeDofs;
  ss << std::endl;
  ss << "Element DoF types:";
  for (uint16 i = 0; i < Element::getNumberOfDofs(); i++) {
    ss << " " << Dof::dofTypeLabels[Element::getDofType(i)];
  }
  ss << " Number of element DoFs: " << numberOfElementDofs;
  ss << std::endl;

  echolog(ss.str().c_str());
	echolog("DoFs = %d, constrained DoFs = %d, MPC eq. = %d, TOTAL eq. = %d",
                n_dofs, n_constrained_dofs, n_MPC_eqs, n_solve_dofs + n_MPC_eqs);

	if (!material) {
		warning("FEStorage::prepare_for_solution: material isn't defined");
		return false;
	}
	return true;
}


// Kij_add is a function to add a value to a global stiffnes matrix 
// from element stiffnes matricies. As long as we have distinct blocks
// of global matrix for constrained DoFs, MPC's lambdas and DoFs to be
// solver we need to choose in which block (Kcc, Kcs, KssCsT)
// the value should be added.
// NOTE: nodes numbers are started from 1. And DoFs are started from 0.
void FEStorage::Kij_add(int32 nodei, Dof::dofType dofi, int32 nodej, Dof::dofType dofj, double value) {
  Dof* rowDof = getDof(nodei, dofi);
  Dof* colDof = getDof(nodej, dofj);
	uint32 rowEq = rowDof->eqNumber;
	uint32 colEq = colDof->eqNumber;
	if (rowEq > colEq) swap(rowEq, colEq);
	if (rowEq <= n_constrained_dofs) {
		if (colEq <= n_constrained_dofs) {
			Kcc->addValue(rowEq, colEq, value);
    } else {
			Kcs->addValue(rowEq, colEq - n_constrained_dofs, value);
    }
  } else {
		KssCsT->addValue(rowEq - n_constrained_dofs, colEq - n_constrained_dofs, value);
  }
}


// Cij_add is a function to add a coefficient from MPC equation to the global
// matrix.
// eq_num - number of MPC equation,
// nodej, dofj - DoFs for which the coefficient to be set.
void FEStorage::Cij_add(uint32 eq_num, int32 nodej, Dof::dofType dofj, double coef) {
	assert(Cc);
	assert(eq_num > 0 && eq_num <= n_MPC_eqs);
	uint32 dof_col = getDofEqNumber(nodej, dofj);
	if (dof_col <= n_constrained_dofs) {
		Cc->addValue(eq_num, dof_col, coef);
  }	else {
		//Cs
		KssCsT->addValue(dof_col-n_constrained_dofs, n_solve_dofs+eq_num, coef);
  }
}

// Fi_add is a function to add the value to a rhs of the global system of linear
// equations. 
// NOTE: nodes index starts from 1. DoF index starts with 0.
void FEStorage::Fi_add(int32 nodei, Dof::dofType dofi, double value) {
	assert(vec_F);
	uint32 row = getDofEqNumber(nodei, dofi);
	assert(row <= n_dofs);
	vec_F[row-1] += value;
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
	assert(vec_F);
	memset(vec_F, 0, sizeof(double)*(n_dofs));
	memset(vec_dq_dlambda, 0, sizeof(double)*(n_dofs+n_MPC_eqs));
  if (n_MPC_eqs > 0) {
    memset(vec_b, 0, sizeof(double)*(n_MPC_eqs));
  }
}

// get_solve_mat:
// just getter
math::SparseSymmetricMatrix& FEStorage::get_solve_mat() {
	assert(KssCsT);
	return *KssCsT;
}

// get_solver_rhs:
// prepare and pass back a rhs of the global eq. system.
// As we exclude constrained DoFs from a eq. system we need to modify rhs of the
// system.
double* FEStorage::get_solve_rhs () {
	assert(vec_rhs);
	double *KcsTdqc = new double[n_solve_dofs];
	Kcs->transpose_mult_vec(vec_dqc,KcsTdqc);
	for (uint32 i=0; i < n_solve_dofs; i++) {
		vec_rhs_dof[i] = vec_Fs[i] - KcsTdqc[i];//Kcs->transpose_mult_vec_i(vec_dqc,i+1);//TODO: тут может быть ошибка
  }
	for (uint32 i=0; i < n_MPC_eqs; i++) {
		vec_rhs_mpc[i] = vec_b[i] - Cc->mult_vec_i(vec_dqc,i+1);
  }
	delete[] KcsTdqc;
	return vec_rhs;
}


uint16 FEStorage::add_post_proc (PostProcessor *pp) {
	assert(pp);
	uint16 num = static_cast<uint16> (this->post_procs.size()+1);
	pp->nPost_proc = num;
	post_procs.push_back(pp);
	return num;
}

// clearMesh:
// delete element and node tables and clear boundary conditions
void FEStorage::clearMesh () {
	status = ST_INIT; //TODO продумать
	deleteElements();
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
}

// deleteElements:
// delete element table
void FEStorage::deleteElements() {
  for (uint32 i = 0; i < n_elements; i++) {
    delete elements[i];
  }
  elements.clear();
  n_elements = 0;
}

//nodes_reassign(_nn)
void FEStorage::nodes_reassign(uint32 _nn)
{
	nodes.clear();
	n_nodes = _nn;
  //Node() fires Vec<3> constructor, thus Node coordinates are (0,0,0) by default
  //TODO: try-catch of memory overflow
	nodes.assign(_nn, Node());
}

//elements_reassign(_en)
void FEStorage::elements_reassign(uint32 _en) {
  deleteElements();
	n_elements = _en;
  //TODO: catch if not enough memory
  elements.reserve(_en);
  ElementFactory::createElements (elType, n_elements, elements); 
  Element::storage = this;
  for (uint32 i = 0; i < _en; i++) {
    //access elNum protected values as friend
    elements[i]->elNum = i+1;
  }
}


// get_q_e(el, ptr) функция возвращает вектор узловых степеней свободы элемента,
// вызывающая сторона должна предоставить массив ptr размерностью Element::n_nodes()*Node::n_dofs() + Element::n_dofs()
// el начинается с 1
//void FEStorage::get_q_e(uint32 el, double* ptr)
//{
//	assert(el <= n_elements);
//	assert(vec_q_lambda);
//	for (uint16 i=0; i<Element::n_nodes(); i++)
//		for (uint16 j=0; j<Node::n_dofs(); j++)
//			ptr[i*Node::n_dofs()+j] = vec_q_lambda[getDofEqNumber(elements[el-1]->getNodeNumber(i), j)-1];
//	for (uint16 i=0; i < Element::n_dofs(); i++)
//		ptr[Element::n_nodes()*Node::n_dofs()+i] = vec_q_lambda[getDofEqNumber(-(int32)el, i)-1];
//}
//// get_dq_e см. выше
//void FEStorage::get_dq_e(uint32 el, double* ptr)
//{
//	assert(el <= n_elements);
//	assert(vec_dq_dlambda);
//	for (uint16 i=0; i<Element::n_nodes(); i++)
//		for (uint16 j=0; j<Node::n_dofs(); j++)
//			ptr[i*Node::n_dofs()+j] = vec_dq_dlambda[getDofEqNumber(elements[el-1]->getNodeNumber(i), j)-1];
//	for (uint16 i=0; i < Element::n_dofs(); i++)
//		ptr[Element::n_nodes()*Node::n_dofs()+i] = vec_dq_dlambda[getDofEqNumber(-(int32)el, i)-1];
//}

double FEStorage::getDofSolution (int32 node, Dof::dofType dof) {
  assert(vec_q_lambda);
  return vec_q_lambda[getDofEqNumber(node, dof) - 1];
}

double FEStorage::getDofSolutionDelta (int32 node, Dof::dofType dof) {
  assert(vec_dq_dlambda);
  return vec_dq_dlambda[getDofEqNumber(node, dof) - 1];
}

bool FEStorage::isDofUsed (int32 node, Dof::dofType dof) {
  Dof* ptr;
  if (node < 0) {
    if (!Element::isDofUsed(dof)) {
      return false;
    }
    ptr = elementDofs[(node-1)*Element::getNumberOfDofs() +Element::getDofIndex(dof)];
    if (ptr == NULL) {
      return false;
    }
  } else {
    if (!Node::isDofUsed(dof)) {
      return false;
    }
    ptr = nodeDofs[(node-1)*Node::getNumberOfDofs() + Node::getDofIndex(dof)];
    if (ptr == NULL) {
      return false;
    }
  }
  return true;
}

// n starts from 1
// prt array always has 3 elements
void FEStorage::getNodePosition(uint32 n, double* ptr, bool deformed)
{
	assert(n > 0 && n <= n_nodes);
	for (uint16 i = 0; i < 3; i++) {
		ptr[i] = nodes[n-1].pos[i];
  }
	if (deformed) {
    if (isDofUsed(n, Dof::UX)) {
      ptr[0] += getDofSolution(n, Dof::UX);
    }

    if (isDofUsed(n, Dof::UY)) {
      ptr[1] += getDofSolution(n, Dof::UY);
    }

    if (isDofUsed(n, Dof::UZ)) {
      ptr[2] += getDofSolution(n, Dof::UZ);
    }
	}
}

double FEStorage::getReaction(int32 n, Dof::dofType dof) {
  if (isDofUsed(n, dof)) {
    uint32 eq_num = getDofEqNumber(n,dof);
    assert(vec_reactions);
    //assert(eq_num > 0 && eq_num <= n_constrained_dofs);
    if (eq_num > 0 && eq_num <= n_constrained_dofs) {
      return vec_reactions[eq_num-1];
    } else {
      return 0.0;
    }
  } else {
    return 0.0;
  }
}


void FEStorage::pre_first() {
	//включаем тренировку матриц
	KssCsT->startTraining();
	if (n_MPC_eqs) Cc->startTraining();
	Kcc->startTraining();
	Kcs->startTraining();
}


void FEStorage::post_first() {
	//выключаем тренировку матриц
	if (n_MPC_eqs) {
		Cc->stopTraining();
	}
	KssCsT->stopTraining();
	
	Kcc->stopTraining();
	Kcs->stopTraining();

}


void FEStorage::process_solution()
{
	// из вектора решения вытаскиваем все необходимое
	//складываем решение

	for (uint32 i=0; i < n_dofs; i++) {
		vec_q[i] += vec_dq[i];
  }
	for (uint32 i=n_dofs; i < n_dofs+n_MPC_eqs; i++) {
		vec_q_lambda[i] = vec_dq_dlambda[i];
  }
	
	for (uint32 i = 0; i < n_constrained_dofs; i++) {
		vec_reactions[i] = Kcs->mult_vec_i(vec_dqs,i+1) + Kcc->mult_vec_i(vec_dqc,i+1) - vec_Fc[i];
		if (n_MPC_eqs && Cc->getNumberOfValues()) {
			vec_reactions[i] += Cc->transpose_mult_vec_i(vec_dlambda,i+1);
    }
	}
}


void FEStorage::apply_BCs (uint16 curLoadstep, uint16 curIteration, double d_par, double cum_par) {
  // terms of MPCs: Cc Cs matricies and vec_b vector
  // because of non-linear MPC we need to update
  // MPC coefficients every step
  //
  for (size_t i = 0; i < mpcCollections.size(); i++) {
    mpcCollections[i]->update();
    //mpcCollections[i]->printEquations(std::cout);
  }

  uint32 eq_num = 1;
  list<Mpc*>::iterator mpc = list_bc_MPC.begin();
  while (mpc != list_bc_MPC.end()) {
    vec_b[eq_num-1] = (*mpc)->b;
    list<MpcTerm>::iterator token = (*mpc)->eq.begin();
    while (token != (*mpc)->eq.end()) {
      Cij_add(eq_num, token->node, token->node_dof, token->coef);
      token++;
    }
    eq_num++;
    mpc++;
  }
	
	// fill nodal forces
	list<BC_dof_force>::iterator bc_force = list_bc_dof_force.begin();
	while (bc_force != list_bc_dof_force.end())	{
		Fi_add(bc_force->node, bc_force->node_dof, bc_force->value*cum_par);
		bc_force++;
	}

	// fill nodal displacements (kinematic constraints)
	list<BC_dof_constraint>::iterator bc_dof = list_bc_dof_constraint.begin();
	while (bc_dof != list_bc_dof_constraint.end()) {
		uint32 eq_num = getDofEqNumber(bc_dof->node, bc_dof->node_dof);//TODO: DEBUG stuff
		vec_dq[eq_num-1] = bc_dof->value*d_par;
		bc_dof++;
	}
}


// read Ansys Mechanical APDL *.cdb file. Nodes, Elements, Displacement BC and MPC (Constraint equations) is supported
// read_ans_data repcales storage's mesh.
bool read_ans_data(const char *filename, FEStorage *storage)
{
	uint32 n_number, en;
	ifstream file(filename);
	if (!file)
	{
		warning("read_ans_data: Can't open input file `%s`",filename);
		return false;
	}
	storage->clearMesh();
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
        error("read_ans_data: NBLOCK: maximum node number is %d, but \
            number of nodes is %s. Note that nla3d needs compressed numbering \
            for nodes and elements", max_n_number, n_number );
      }
			storage->nodes_reassign(n_number);
			file.getline(buf, 1024);
      string buf_str(buf);
      // we need to take a format of columns "3i9"
      size_t start = buf_str.find("i")+1;
      size_t  stop = buf_str.find(",");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str());

      //TODO: as we need numbering from 1/0 to n_number, here we can check that number of nodes and theirs id are in a row
			for (uint32 i=1; i<=n_number; i++)
			{
				file.getline(buf, 1024);
				size_t len=strlen(buf);
				for (uint16 j=0; j<3;j++)
					if (len>=3*frmt+20*(j+1))
            //note that last column in NBLOCK table could be avoided if Z=0.0
            //but storage->nodes_reassign(n_number) initialize the node table with (0,0,0)
						storage->getNode(i).pos[j] = atof(string((char*) (buf+3*frmt+20*j),20).c_str());
			}
		}//NBLOCK
		else if (vec[0].find("EBLOCK")!=vec[0].npos)
		{  
      //EBLOCK,19,SOLID,      7024,      7024
      //(19i9)
			en = atoi(vec[6].c_str());
      if (en != atoi(vec[8].c_str())) {
        error("read_ans_data: EBLOCK: maximum element number is %d, but \
            number of element is different. Note that nla3d needs compressed numbering \
            for nodes and elements", en);
      }
			storage->elements_reassign(en);
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
//TODO: not a good thing.. It seems that getline keeps windows line ending
//#ifdef linux
//        len--;
//#endif
        //TODO: check that n_nodes and provided number of nodes are the same
        if (len != 11*frmt+frmt*Element::n_nodes())
          warning("read_ans_data: in EBLOCK for element %d the number of nodes provided is not equal to %d", i, Element::n_nodes());
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
        storage->add_bounds(bnd);
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
			storage->add_bounds(mpc);
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
            for (size_t i = rangesVec[all-1] + 1; i < -rangesVec[all]+1; i++) {
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
	storage->setStatus(ST_LOADED);
	return true;
}



uint32 FEStorage::getDofEqNumber(int32 node, Dof::dofType dof) {
  return getDof(node, dof)->eqNumber;
}

Dof* FEStorage::getElementDof (int32 el, Dof::dofType dof) {
  Dof* ptr = elementDofs[(el-1)*Element::getNumberOfDofs() + Element::getDofIndex(dof)];
  assert(ptr != NULL);
  return ptr;
}

Dof* FEStorage::getNodeDof (int32 node, Dof::dofType dof) {
  Dof* ptr = nodeDofs[(node-1)*Node::getNumberOfDofs() + Node::getDofIndex(dof)];
  assert(ptr != NULL);
  return ptr;
}

Dof* FEStorage::getDof (int32 node, Dof::dofType dof) {
  if (node < 0) {
    return getElementDof(-node, dof);
  }
  return getNodeDof(node, dof);
}

// element_nodes(el, node_ptr), вызывающая сторона должна предоставить массив 
// указателей Node* на >= Element::nNodes() элементов
// el начинается с 1
void FEStorage::element_nodes(uint32 el, Node** node_ptr)
{
	assert(el <= n_elements);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		node_ptr[i] = & nodes[elements[el-1]->getNodeNumber(i)-1];
}

void FEStorage::addFEComponent (FEComponent* comp) {
	assert(comp);
	feComponents.push_back(comp);
}

FEComponent* FEStorage::getFEComponent(size_t i) {
  assert(i < feComponents.size());
  return feComponents[i];
}

FEComponent* FEStorage::getFEComponent(string name) {
  for (size_t i = 0; i < feComponents.size(); i++) {
    if (name.compare(feComponents[i]->name) == 0) {
      return feComponents[i];
    }
  }
  error("FEStorage::FEComponent: Can't find a component with name %s", name.c_str());
}

void FEStorage::listFEComponents (std::ostream& out) {
  for (size_t i = 0; i < feComponents.size(); i++) {
    feComponents[i]->print(out);
    out << std::endl;
  }
}

void FEStorage::registerNodeDof(int32 node, Dof::dofType dof) {
  if (nodeDofs.size() == 0) {
    nodeDofs.assign(getNumberOfNodes()*Node::getNumberOfDofs(), NULL);
  }
  uint32 ind = (node-1)*Node::getNumberOfDofs() + Node::getDofIndex(dof);
  if (nodeDofs[ind] == NULL) {
    nodeDofs[ind] = new Dof;
    numberOfNodeDofs++;
  }
}

void FEStorage::registerElementDof(int32 el, Dof::dofType dof) {
  if (elementDofs.size() == 0) {
    elementDofs.assign(getNumberOfElements()*Element::getNumberOfDofs(), NULL);
  }
  uint32 ind = (el-1)*Element::getNumberOfDofs() + Element::getDofIndex(dof);
  if (elementDofs[ind] == NULL) {
    elementDofs[ind] = new Dof;
    numberOfElementDofs++;
  }
}


} // namespace nla3d 

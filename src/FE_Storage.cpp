#include "FE_Storage.h"

FE_Storage::FE_Storage()  {
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

	dof_array = NULL;
	
	status = ST_INIT;
  elType = ElementFactory::NOT_DEFINED;
};

FE_Storage::~FE_Storage () {
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

  clearMesh();
}

bool FE_Storage::prepare_for_solution () {
	if (status < ST_LOADED) {
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


//функция возвращает ссылку на значение, которое хранится в глобальной матрице жосткости по адресу 
// строка: dofi-тая степень свободы nodei-того узла 
// столбец: dofj-тая степень свободы nodej-того узла 
// узлы нумеруются с индекса 1, dofi=[0;Node::nDOF()-1];
void FE_Storage::Kij_add(int32 nodei, uint16 dofi, int32 nodej, uint16 dofj, double value) {
	//assert(matK);
	uint32 eq_row = get_dof_eq_num(nodei, dofi);
	uint32 eq_col = get_dof_eq_num(nodej, dofj);
	if (eq_row > eq_col) swap(eq_row, eq_col);
	if (eq_row <= n_constrained_dofs) {
		if (eq_col <= n_constrained_dofs) {
			Kcc->add_value(eq_row, eq_col, value);
    } else {
			Kcs->add_value(eq_row, eq_col - n_constrained_dofs, value);
    }
  } else {
		KssCsT->add_value(eq_row - n_constrained_dofs, eq_col - n_constrained_dofs, value);
  }
}

void FE_Storage::Cij_add(uint32 eq_num, int32 nodej, uint32 dofj, double value) {
	assert(Cc);
	assert(eq_num > 0 && eq_num <= n_MPC_eqs);
	uint32 dof_col = get_dof_eq_num(nodej, dofj);
	if (dof_col <= n_constrained_dofs) {
		Cc->add_value(eq_num, dof_col, value);
  }	else {
		//Cs
		KssCsT->add_value(dof_col-n_constrained_dofs, n_solve_dofs+eq_num, value);
  }
}

//getFi(..) возвращает ссылку на ячейку глобального вектора сил
//по адресу dofi-тая степень свободы nodei-того узла
// узлы нумеруются с индекса 1, dofi=[0;Node::nDOF()-1];
// TODO: CHECK
void FE_Storage::Fi_add(int32 nodei, uint16 dofi, double value) {
	assert(vec_F);
	uint32 row = get_dof_eq_num(nodei, dofi);
	assert(row <= n_dofs);
	vec_F[row-1] += value;
}

void FE_Storage::zeroK() {
	assert(KssCsT);
	KssCsT->zero_block(n_solve_dofs);// TODO: сделать эту функцию //DEBUG
	//KssCsT->zero();//DEBUG
	Kcc->zero();
	Kcs->zero();
}

void FE_Storage::zeroF() {
	assert(vec_F);
	memset(vec_F, 0, sizeof(double)*n_dofs);
	memset(vec_dq_dlambda, 0, sizeof(double)*(n_dofs+n_MPC_eqs));//TODO: test
}


Sparse_Matrix_SUR& FE_Storage::get_solve_mat() {
	assert(KssCsT);
	return *KssCsT;
}


double* FE_Storage::get_solve_rhs () {
	assert(vec_rhs);
	double *KcsTdqc = new double[n_solve_dofs];
	Kcs->transpose_mult_vec(vec_dqc,KcsTdqc);
	for (uint32 i=0; i < n_solve_dofs; i++)
		vec_rhs_dof[i] = vec_Fs[i] - KcsTdqc[i];//Kcs->transpose_mult_vec_i(vec_dqc,i+1);//TODO: тут может быть ошибка
	for (uint32 i=0; i < n_MPC_eqs; i++)
		vec_rhs_mpc[i] = vec_b[i] - Cc->mult_vec_i(vec_dqc,i+1);
	delete[] KcsTdqc;
	return vec_rhs;
}


uint16 FE_Storage::add_post_proc (Post_proc *pp) {
	assert(pp);
	uint16 num = this->post_procs.size()+1;
	pp->nPost_proc = num;
	post_procs.push_back(pp);
	return num;
}

//clearMesh() функция очищает таблицу узлов, элементов, ГУ
void FE_Storage::clearMesh () {
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


void FE_Storage::deleteElements() {
  for (uint32 i = 0; i < n_elements; i++) {
    delete elements[i];
  }
  elements.clear();
  n_elements = 0;
}

//nodes_reassign(_nn)
void FE_Storage::nodes_reassign(uint32 _nn)
{
	nodes.clear();
	n_nodes = _nn;
  //Node() fires Vec<3> constructor, thus Node coordinates are (0,0,0) by default
  //TODO: try-catch of memory overflow
	nodes.assign(_nn, Node()); // TODO: тут бы catch на возможность выделения памяти
}

//elements_reassign(_en)
void FE_Storage::elements_reassign(uint32 _en)
{
  deleteElements();
	n_elements = _en;
  //TODO: use vector : elements.reserve
  elements.reserve(_en);
  //elements = new Element*[n_elements];
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
void FE_Storage::get_q_e(uint32 el, double* ptr)
{
	assert(el <= n_elements);
	assert(vec_q_lambda);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		for (uint16 j=0; j<Node::n_dofs(); j++)
			ptr[i*Node::n_dofs()+j] = vec_q_lambda[get_dof_eq_num(elements[el-1]->node_num(i), j)-1];
	for (uint16 i=0; i < Element::n_dofs(); i++)
		ptr[Element::n_nodes()*Node::n_dofs()+i] = vec_q_lambda[get_dof_eq_num(-(int32)el, i)-1];
}
// get_dq_e см. выше
void FE_Storage::get_dq_e(uint32 el, double* ptr)
{
	assert(el <= n_elements);
	assert(vec_dq_dlambda);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		for (uint16 j=0; j<Node::n_dofs(); j++)
			ptr[i*Node::n_dofs()+j] = vec_dq_dlambda[get_dof_eq_num(elements[el-1]->node_num(i), j)-1];
	for (uint16 i=0; i < Element::n_dofs(); i++)
		ptr[Element::n_nodes()*Node::n_dofs()+i] = vec_dq_dlambda[get_dof_eq_num(-(int32)el, i)-1];
}

//get_q_n(n, ptr)
// n начинается с 1
void FE_Storage::get_q_n(uint32 n, double* ptr)
{
	assert(n > 0 && n <= n_nodes);
	assert(vec_q_lambda);
	for (uint16 j=0; j<Node::n_dofs(); j++)
		ptr[j] = vec_q[get_dof_eq_num(n, j)-1];
}

//void get_node_pos(uint32 n, double* ptr, bool def = false) double массим на 3 элемента!
//n с 1
void FE_Storage::get_node_pos(uint32 n, double* ptr, bool def)
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

double FE_Storage::get_reaction_force(int32 n, uint16 dof) {
	uint32 eq_num = get_dof_eq_num(n,dof);
	assert(vec_reactions);
	assert(eq_num > 0 && eq_num <= n_constrained_dofs);
	return vec_reactions[eq_num-1];
}


void FE_Storage::pre_first() {
	//включаем тренировку матриц
	KssCsT->start_training();
	if (n_MPC_eqs) Cc->start_training();
	Kcc->start_training();
	Kcs->start_training();
}


void FE_Storage::post_first() {
	//выключаем тренировку матриц
	if (n_MPC_eqs) 
	{
		Cc->stop_training();
	}
	KssCsT->stop_training();
	
	Kcc->stop_training();
	Kcs->stop_training();

}


void FE_Storage::process_solution()
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


void FE_Storage::apply_BCs (uint16 curLoadstep, uint16 curIteration, double d_par, double cum_par) {
	if (curLoadstep == 1 && curIteration == 1) {
		//заполним Cc Cs и vec_b
		uint32 eq_num = 1;
		list<BC_MPC>::iterator mpc = list_bc_MPC.begin();
		while (mpc != list_bc_MPC.end()) {
			vec_b[eq_num-1] = mpc->b;
			list<MPC_token>::iterator token = mpc->eq.begin();
			while (token != mpc->eq.end()) {
				Cij_add(eq_num, token->node, token->node_dof, token->coef);
				token++;
			}
			eq_num++;
			mpc++;
		}
	}
	
	//заполним вектор узловых сил
	list<BC_dof_force>::iterator bc_force = list_bc_dof_force.begin();
	while (bc_force != list_bc_dof_force.end())	{
		Fi_add(bc_force->node, bc_force->node_dof, bc_force->value*cum_par);
		bc_force++;
	}

	//заполним вектор заданных перемещений
	list<BC_dof_constraint>::iterator bc_dof = list_bc_dof_constraint.begin();
	while (bc_dof != list_bc_dof_constraint.end()) {
		uint32 eq_num = get_dof_eq_num(bc_dof->node, bc_dof->node_dof);//TODO: DEBUG stuff
		vec_dq[eq_num-1] = bc_dof->value*d_par;
		bc_dof++;
	}
}


// read Ansys Mechanical APDL *.cdb file. Nodes, Elements, Displacement BC and MPC (Constraint equations) is supported
// read_ans_data repcales storage's mesh.
bool read_ans_data(const char *filename, FE_Storage *storage)
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
			n_number= atoi(vec[6].c_str());
      if (max_n_number != n_number) {
        warning("read_ans_data: NBLOCK: maximum node number is %d, but number of nodes is %s. Note that nla3d needs compressed numbering for nodes and elements", max_n_number, n_number );
      }
			storage->nodes_reassign(n_number);
			file.getline(buf, 1024);
      string buf_str(buf);
      // we need to take a format of columns "3i9"
      uint16 start = buf_str.find("i")+1;
      uint16 stop = buf_str.find(",");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str());

      //TODO: as we need numbering from 1/0 to n_number, here we can check that number of nodes and theirs id are in a row
			for (uint32 i=1; i<=n_number; i++)
			{
				file.getline(buf, 1024);
				uint16 len=strlen(buf);
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
			en=atoi(vec[6].c_str());
			storage->elements_reassign(en);
			file.getline(buf, 1024);
      // we need to take a format of columns "3i9"
      // in Ansys 12 here is 8 symbols per number (19i8), but in ansys 15 (19i9) is used. 
      string buf_str(buf);
      uint16 start = buf_str.find("i")+1;
      uint16 stop = buf_str.find(")");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str()); 
			for (uint32 i=1; i<=en; i++)
			{
				file.getline(buf, 1024);
				uint16 len=strlen(buf);
        //TODO: check that n_nodes and provided number of nodes are the same
        if (len != 11*frmt+frmt*Element::n_nodes())
          warning("read_ans_data: in EBLOCK for element %d the number of nodes provided is not equal to %d", i, Element::n_nodes());
				for (uint16 j=0; j<Element::n_nodes();j++)
					if (len>=11*frmt+frmt*(j+1))
						storage->getElement(i).node_num(j) = atoi(string((char*) (buf+11*frmt+frmt*j),frmt).c_str());
			}
		}//EBLOCK
		else if (vec[0].find('D')!=vec[0].npos && vec[0].length()==1)
		{
				BC_dof_constraint bnd;
				bnd.node = atoi(vec[2].c_str());
				bnd.value = atof(vec[6].c_str());
        bnd.node_dof = str2dof(vec[4]);
        storage->add_bounds(bnd);
		}//D
    else if (vec[0].compare("CE") == 0)
    {
      //How MPC looks like this in cdb file:
      //CE,R5.0,DEFI,       2,       1,  0.00000000    
      //CE,R5.0,NODE,      1700,UX  ,  1.00000000    ,      1700,UZ  ,  1.00000000  
      BC_MPC mpc;
      mpc.b = atoi(vec[10].c_str()); //rhs of MPC equation
      uint16 n_terms = atoi(vec[6].c_str()); //number of terms in equation
      //debug("MPC link: %d terms, b = %f", n_terms, mpc.b);
      while (n_terms > 0)
      {
        file.getline(buf, 1024);
        vector<string> vec = read_tokens(buf);
        uint16 place = 6;
        for (int i=0; i < max(n_terms, 2); i++) 
        {
          uint32 node = atoi(vec[place].c_str());
          uint16 dof = str2dof(vec[place+2]);
          double coef = atof(vec[place+4].c_str());
          //debug("%d term: node %d, dof %d, coef = %f", i, node, dof, coef);
          mpc.eq.push_back(MPC_token(node,dof,coef));
          place += 6;
          n_terms--;
        }
      }
			storage->add_bounds(mpc);
    }//CE (MPC)
    //TODO: add FX FY FZ
	}
	file.close();
	storage->setStatus(ST_LOADED);
	return true;
}


//I'd like to make this functions inline

uint32 FE_Storage::get_dof_num(int32 node, uint16 dof) {
	// возвращает число от 1 до n_dofs
	uint32 res = (node < 0)?((-node-1)*Element::n_dofs()+dof+1):(n_elements*Element::n_dofs()+(node-1)*Node::n_dofs()+dof+1);
	return res; 
}

uint32 FE_Storage::get_dof_eq_num(int32 node, uint16 dof) {
	assert(dof_array);
	return dof_array[get_dof_num(node,dof)-1].eq_number;
}

// element_nodes(el, node_ptr), вызывающая сторона должна предоставить массив 
// указателей Node* на >= Element::nNodes() элементов
// el начинается с 1
void FE_Storage::element_nodes(uint32 el, Node** node_ptr)
{
	assert(el <= n_elements);
	for (uint16 i=0; i<Element::n_nodes(); i++)
		node_ptr[i] = & nodes[elements[el-1]->node_num(i)-1]; //TODO: CHECK
}

bool FE_Storage::is_dof_constrained(int32 node, uint16 dof) {
	assert(dof_array);
	return dof_array[get_dof_num(node,dof)-1].is_constrained;
}



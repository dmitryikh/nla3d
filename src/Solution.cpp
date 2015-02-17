#include "Solution.h"



bool Solution::setInit ()
{
	if (status != SOL_SOLVING)
	{
		qLoadstep = 10;
		qIterat = 20;
		Criteria = 1e-3;
		solver_type = SAE_CHOLESKY;
		return true;
	}
	warning("Solution::setInit: can't change solution parameters while solving");
	return false;
}
bool Solution::setqIterat (uint16 iter)
{
	if (status != SOL_SOLVING)
	{
		if (iter < 1)
		{
			warning("Solution::setqIter: qIterat can't be less than 1 (input %d)",iter);
			return false;
		}
		qIterat = iter;
		return true;
	}
	warning("Solution::setqIter: can't change solution parameters while solving");
	return false;
}
bool Solution::setqLoadstep (uint16 ls)
{
	if (status != SOL_SOLVING)
	{
		if (ls < 1)
		{
			warning("Solution::setqLoadstep: qLoadstep can't be less than 1 (input %d)",ls);
			return false;
		}
		qLoadstep = ls;
		return true;
	}
	warning("Solution::setqLoadstep: can't change solution parameters while solving");
	return false;
}
bool Solution::setCriteria (double cr)
{
	if (status != SOL_SOLVING)
	{
		if (cr < 0.0f)
		{
			warning("Solution::setCriteria: Criteria must be positive (input %f)",cr);
			return false;
		}
		Criteria = cr;
		return true;
	}
	warning("Solution::setCriteria: can't change solution parameters while solving");
	return false;
}
void Solution::main_process ()
{
	echolog("Start the solution process");
	assert(storage);
	Timer sol(true);
	Timer overall(true);
	if (!storage->prepare_for_solution())
	{
		warning("Solution::main_process: solution is interrupted");
		status = SOL_WAIT;
		return;
	}

	for (uint16 i=0; i < storage->getNumPostProc(); i++)
		storage->getPostProc(i).pre(qLoadstep);

	for (uint32 el = 1; el <= storage->getNumElement(); el++)
				storage->getElement(el).pre(el, storage);
	echolog("Preporation time: %4.2f sec.", sol.stop());
	curCriteria = 0.0f;
	double dF_par;
	double cumF_par = 0.0;
	uint16 cum_iterations = 0;
	for (curLoadstep = 1; curLoadstep <= qLoadstep; curLoadstep++ )
	{
		dF_par = 1.0/qLoadstep;
		cumF_par += dF_par;
		
		curIterat= 1;
		for (; curIterat <= qIterat; curIterat++ )
		{
			cum_iterations++;
			echolog("------LS = %d of %d, IT = %d (%d)------", curLoadstep, qLoadstep, curIterat, cum_iterations);
			if (stopit)
			{
				curLoadstep--;
				echolog("Solution::main_process: solution process was stopped by user.\n Last solved loadstep = %d",curLoadstep);
				status = SOL_WAIT;
				return;
			}
			if (curIterat > 1) dF_par = 0.0;
			if (cum_iterations == 1)
				storage->pre_first();
			else
			{
				// очистка matK и vecF
				storage->zeroK();
				storage->zeroF();
			}
			if (cum_iterations == 1)
			{
				sol.start();
				echolog("Start formulation for elements ( %d )", storage->getNumElement());
			}
			for (uint32 el = 1; el <= storage->getNumElement(); el++)
				storage->getElement(el).build(el, storage);
			if (cum_iterations == 1)
				echolog("Formulation time: %4.2f sec.", sol.stop());
			//учет "узловых" Г.У. пока только кинематические
			if (cum_iterations == 1)
				sol.start();
			storage->apply_BCs(curLoadstep, curIterat, dF_par, cumF_par);
			if (cum_iterations == 1)
			{
				echolog("BC applying time: %4.2f sec.", sol.stop());
				sol.start();
				storage->post_first();
				echolog("Matrices training time: %4.2f sec.", sol.stop());
			}
			//решение системы
			if (!solve_wrap())
			{
				warning ("Solution::main_process: solver failed. exiting");
				status = SOL_WAIT;
				return;
			}
			if (cum_iterations == 1)
				sol.start();
			storage->process_solution();
			if (cum_iterations == 1)
				echolog("After solution process time: %4.2f sec.", sol.stop());

			if (cum_iterations == 1)
				sol.start();
			//вычисление деформаций и напряжений в элементах
			for (uint32 el = 1; el <= storage->getNumElement(); el++)
				storage->getElement(el).update(el, storage);

			if (cum_iterations == 1)
				echolog("Elements updating time: %4.2f sec.", sol.stop());

			curCriteria = 0.0;
			for (uint32 i=0; i < storage->get_n_solve_dofs(); i++)
			{
				curCriteria += fabs(storage->get_vec_dqs()[i]);
			}
			curCriteria /= storage->get_n_solve_dofs();

			echolog("Criteria (%5.1f%%)", curCriteria/Criteria*100);
			if (curCriteria < Criteria) 
			{
				curIterat++;
				break;
			}
			
		}//iterations
		echolog("LS %d completed with %d IT", curLoadstep, curIterat-1);
		if (curLoadstep == 1)
			sol.start();
		for (uint16 i=0; i < storage->getNumPostProc(); i++)
			storage->getPostProc(i).process(curLoadstep, qLoadstep);
		if (curLoadstep == 1)
			echolog("PostProc's processing time: %4.2f sec.", sol.stop());
		echolog("\n\n");
	}//LS steps
	storage->setStatus(ST_SOLVED);
	status = SOL_WAIT;
	echolog("****SOLVED****  (Total solution time: %4.2f sec.)", overall.stop());
	for (uint16 i=0; i < storage->getNumPostProc(); i++)
			storage->getPostProc(i).post(curLoadstep, qLoadstep);
}
uint16 Solution::solve_wrap ()
{
	/*uint16 res;
	solver_type = SAE_CHOLESKY;
	if (solver_type == SAE_CHOLESKY)
	{
		res = solveSAE_Cholesky();
		if (!res)
		{
			warning("Switch solver to Bunch-Kaufman factorization..");
			solver_type = SAE_BUNCH;
			res = solveSAE_Bunch();
			if (!res)
				error("Solution process is terminated");
		}
	}
	else if(solver_type == SAE_BUNCH)
	{
		res = solveSAE_Bunch();
		if (!res)
		error("Solution process is terminated");
	}
	else
		warning("Unknown solver!");*/
	return solveSAE_DSS();
}
// SAE_PS - solver for positive def. symm. matrix
uint16 Solution::solveSAE_Cholesky()
{
	////Внимание в процессе решения vecF изменяется
	////и уже не будут содежрать информации о нагрузке
	//int n = storage->getNumDofs();	//TODO: перевод в int недопустим при больших размерностях
	//int kd = storage->get_nBand() - 1; //^
	//int ldab = kd+1;	//^
	//int info;
	//int nrhs = 1;
	//int ldb = n;	//^

	//Mat_Band_cm matrix(n, storage->get_nBand());
	//matrix = storage->getMatK();
	//dpbtrf("U",&n,&kd,matrix.Ptr(),&ldab,&info);
	//if (info < 0) error("Solution::solve: error in MKL solver (errno %d)", info);
	//if (info > 0)
	//{
	//	warning("Solution::solve: MKL POSITIVE solver: the leading minor of order %d isn`t positive defined",info);
	//		return 0;
	//}
	//dpbtrs("U",&n,&kd,&nrhs,matrix.Ptr(),&ldab,storage->getVecF().Ptr(),&ldb,&info);
	//if (info != 0) error("Solution::solve: error in MKL solver (errno %d)", info);
	//storage->getVecdQ() = storage->getVecF();// TODO CHECK
	return 1;
}


uint16 Solution::solveSAE_Bunch ()
{
	////Внимание в процессе решения матрица matK и vecF изменяются
	////и уже не будут содежрать информации о матрице жест. и нагрузке
	//int n = storage->getNumDofs();	//TODO: перевод в int недопустим при больших размерностях
	//int info;
	//int nrhs = 1;
	//int ldb = n;	//^

	//Vec_Long triag_m(n*(n+1)/2);
	////triag_m.zero();
	//int *ipv= new int[n];
	//storage->getMatK().toUTriangular(triag_m.Ptr());
	//dsptrf("U", &n, triag_m.Ptr(), ipv, &info);
	//if (info < 0) error("Solution::solve: error in MKL Bunch-Kaufman solver (errno %d)", info);
	//if (info > 0)
	//{
	//	warning("Solution::solve: MKL Bunch-Kaufman solver: the (%d, %d) is 0!",info, info);
	//	return 0;
	//}
	//dsptrs("U", &n, &nrhs,  triag_m.Ptr(), ipv, storage->getVecF().Ptr(), &ldb, &info);
	//if (info != 0) error("Solution::solve: error in MKL Bunch-Kaufman solver (errno %d)", info);
	//delete[] ipv;	
	////storage->vecdQ->operator=(*storage->vecF);
	//storage->getVecdQ() = storage->getVecF();// TODO CHECK
	return 1;
}

uint16 Solution::solveSAE_DSS ()
{
	int n = storage->get_solve_mat().get_n_rows();
	Timer sol_timer;
	int nrhs = 1; // кол-во правых частей
	int mtype = -2; //real symmetric undifinite defined matrix
	int iparm[64];
	int maxfct, mnum, phase, error, msglvl;
	for (uint16 i=0;i<64;i++) iparm[i]=0;
	iparm[0] = 1; //no solver default
	iparm[1] = 2; //fill-in reordering from meris
	iparm[2] = MKL_Get_Max_Threads();
	iparm[3] = 0; //no iterative-direct algorithm
	iparm[4] = 0; //no user fill-in reducing permutation
	iparm[5] = 0; //write solution into x
	iparm[6] = 16; //default logical fortran unit number for output
	iparm[7] = 2; //max numbers of iterative refinement steps
	iparm[9] = 13; //pertrub the pivor elements with 1e-13
	iparm[10] = 1; //use nonsymmetric permutation  and scaling MPS
	iparm[13]=0; //output: number of perturbed pivots
	iparm[17]=-1; //output: number of nonzeros in the factor LU
	iparm[18]=-1; //output: MFLOPS for LU factorization
	iparm[19] = 0; //output: number of CG Iterations

	maxfct = 1; //maximum number of numerical factorizations
	mnum = 1; //wich factorization to use
	msglvl = 0; //dont print statistical information in file
	error = 0; //init error flag
	if (if_solver_first_time)
	{
		sol_timer.start();
		for (uint16 i=0;i<64;i++) pt[i]=0;
		//Reordering and Symbolic Factorization. This step also allocate 
		// all memory that is necessary for the factorization
		phase = 11;
		PARDISO(pt, &maxfct, &mnum, &mtype,&phase,
			&n, storage->get_solve_mat().get_values_array(), (int*) storage->get_solve_mat().get_iofeir_array(), (int*) storage->get_solve_mat().get_columns_array(), 
			NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
		if (error != 0)
		{
			echolog("ERROR during symbolic factorization: %d", error);
			return 0;
		}
		echolog("Number of nonzeros in factors = %d,  Number of factorization MFLOPS = %d", iparm[17],iparm[18]);
		
	}

	// Numerical factorization
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype,&phase,
			&n, storage->get_solve_mat().get_values_array(), (int*) storage->get_solve_mat().get_iofeir_array(), (int*) storage->get_solve_mat().get_columns_array(),
			NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
	if (error != 0)
	{
		echolog("ERROR during numerical factorization: %d", error);
		return 0;
	}

	//Back substitution and iterative refiment
	phase = 33;
	iparm[7] = 2; //max num of iterative refiment
	PARDISO(pt, &maxfct, &mnum, &mtype,&phase,
			&n, storage->get_solve_mat().get_values_array(), (int*) storage->get_solve_mat().get_iofeir_array(), (int*) storage->get_solve_mat().get_columns_array(),
			NULL, &nrhs, iparm, &msglvl, storage->get_solve_rhs(), storage->get_solve_result_vector(), &error);
	if (error != 0)
	{
		echolog("ERROR during solution: %d", error);
		return 0;
	}
	if (if_solver_first_time)
	{
		echolog("Eq. solution time: %f sec.", sol_timer.stop());
		if_solver_first_time = false;
	}
	//Termination and release of memory
	//phase = 1;
	//PARDISO(pt, &maxfct, &mnum, &mtype,&phase,
	//		&n, NULL, (int*) storage->getMatK().get_iofeir_array(), (int*) storage->getMatK().get_columns_array(),
	//		NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);

	//echolog("\nSolve completed ... ");

	return 1;
}
void Solution::bound (double dF_par)
{
	//static const double REAL_BIG_NUMBER = 1.0e10;
	//for(uint32 i=0; i < storage->getNumBound(); i++)
	//{
	//	Bound bnd = storage->getBound(i);
	//	if (bnd.key == D_UX || bnd.key == D_UY || bnd.key == D_UZ)
	//	{
	//		storage->getKij(bnd.node, bnd.key-1,bnd.node, bnd.key-1) *= REAL_BIG_NUMBER;
	//		storage->getFi(bnd.node, bnd.key-1) = storage->getKij(bnd.node, bnd.key-1,bnd.node, bnd.key-1)*bnd.value*dF_par;
	//	}
	//	//TODO: forces dont work now
	//	//else if (storage->bounds[i].key == F_X || storage->bounds[i].key == F_Y || storage->bounds[i].key == F_Z)
	//	//	(*storage->vecF)[row] += storage->bounds[i].value*dF_par;
	//}
}

void Solution::attach (FE_Storage_Interface *st)
{
	assert(st);
	if (status == SOL_SOLVING)
	{
		warning("Solution::attach: can't attach new storage while solving previous one");
		return;
	}
	if (storage) warning("Solution::attach: storage already is attached. Lost the old one");
	storage = st;
}

void Solution::run (bool multi_thread) 
{
	if (status == SOL_SOLVING)
	{
		warning("Solution::run: solution process is already running");
		return;
	}
	status = SOL_SOLVING;
	stopit = false;
	//TODO: кастыль
	if (multi_thread)
		_beginthread(fork,0,this);
	else
		main_process();
}

void Solution::fork (void *ptr)
{
	((Solution*)ptr)->main_process();
}
void Solution::stop ()
{
	stopit = true;
}

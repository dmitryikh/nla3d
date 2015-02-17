#include <conio.h>
#include "sys.h"
#include "math\Mat.h"
#include "math\Mat_Band.h"
#include "FE_Storage.h"
#include "vtk_proc.h"
#include "basic_elements.h"
#include "Solution.h"
#include <mkl.h>
#include "math\Sparse_Matrix.h"
#include <limits>
#include "Reaction_proc.h"

dMat test_func_dmat();
void Mat_test ();
void Sparse_Mat_test();
/*void MKL_test ();
void make_data (FE_Storage &st);
void make_test_data(FE_Storage &st)*/;
void main ()
{
	echolog("---=== WELCOME TO NLA PROGRAM ===---");
	Timer pre_solve(true);
	FE_Storage<MIXED_8N_3D_P0_NL> storage;
	read_ans_data("model.cdb", &storage);
	Material_Comp_Neo_Hookean* mat = new Material_Comp_Neo_Hookean();
	//Material_Hookean* mat = new Material_Hookean();
	mat->Ci(0) = 30;
	mat->Ci(1) = 0.499;
	storage.material = mat;
	Vtk_proc* vtk = new Vtk_proc(&storage);
	Solution sol;
	sol.attach(&storage);
	sol.setqIterat(15);
	sol.setqLoadstep(10);
	double a1 = 130*M_PI/180;
	double a2 = 50*M_PI/180;
	ifstream in("LIST.txt");
	uint32 node;
	double cx,cz;
	double fi_sec = 2*M_PI/3;
	cx = cos(fi_sec+M_PI/2);
	cz = -sin(fi_sec+M_PI/2);
	BC_MPC mpc;
	mpc.b = 0.0;
	while (in)
	{	
		mpc.eq.clear();
		node = 0;
		in >> node;
		if (node)
		{
			mpc.eq.push_back(MPC_token(node,0,cx));
			mpc.eq.push_back(MPC_token(node,2,cz));
			//storage.add_bounds(mpc);
		}
	}
	in.close();
	Reaction_proc proc(&storage);
	list<BC_dof_constraint> &lst = storage.get_dof_const_BC_list();
	list<BC_dof_constraint>::iterator bc_dof = lst.begin();
	while (bc_dof != lst.end())
	{
		if (fabs(bc_dof->value) > 0)
		{
			proc.nodes.push_back(bc_dof->node);
			proc.dofs.push_back(bc_dof->node_dof);
		}
		bc_dof++;
	}
	/*BC_MPC mpc1;
	mpc1.b = 0.0;
	mpc1.eq.push_back(MPC_token(5,0,cos(a1)));
	mpc1.eq.push_back(MPC_token(5,1,sin(a1)));
	BC_MPC mpc2;
	mpc2.b = 0.3;
	mpc2.eq.push_back(MPC_token(3,0,cos(a2)));
	mpc2.eq.push_back(MPC_token(3,1,sin(a2)));*/
	//storage.add_bounds(mpc1);
	//storage.add_bounds(mpc2);
	echolog("Preprocessor time: %f sec.", pre_solve.stop());
	sol.run();
	//Mat_test();
	/*Cmd_Shell *shell = new Cmd_Shell();
	uint32 t1 = tick();
	shell->run();
	delete shell;
	long t2 = clock();
	echolog("Program is closed (Worktime= %10.2fs)",(float) (t2-t1)/1000.0f);*/
	return;
}

void Mat_test ()
{
	Mat<3,3> mat(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f);
	Mat<3,3> mat2 = mat;
	cout << "mat=mat2:" << endl;
	cout << mat << endl;
	cout << "mat2.zero()" << endl;
	mat2.zero();
	cout << mat2 << endl;
	cout << "mat2.Identity()" << endl;
	mat2.Identity();
	mat=mat2;
	cout << mat2 << endl;
	cout << "mat*mat2" << endl;
	cout << mat*mat2 << endl;
	cout << "mat*mat" << endl;
	cout << mat*mat << endl;
	Mat<2,3> mat23(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f);
	cout << "mat23:" << endl << mat23;
	cout << "mat23*mat" << endl << mat23*mat;
	cout << "mat.transpose()" << endl << mat.transpose();
	Vec<3> vec(2.0f, 0.0f, 1.0f);
	cout << "vec:" << endl << vec << endl;
	cout << "mat*vec" << endl <<mat*vec << endl;
	cout << "mat" << endl << mat << endl;
	cout << "mat2" << endl << mat2 << endl;
	cout << "mat+mat" << endl <<mat+mat << endl;
	cout << "mat += mat2" << endl;
	mat += mat2;
	cout << mat << endl;
	//работает

	const Mat<3,3> cmat = mat;
	cout << "cmat:" << endl << cmat << endl;
	cout << "cmat[1][1]: " << cmat[1][1] << endl << endl;

	Mat<3,3> mat33(1.0, 2.0, 3.0, 43.0, 5.0, 6.0, 7.0, 8.0, 20.0);
	Mat<4,4> mat44(1.0, 2.0, 3.0, 4.0, -5.0, 6.0, 7.0, 80.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 20.0);
	Mat<5,5> mat55(1.0, 1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 81.0, 9.0, 12.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 3.0, 18.0, 19.0, 20.0, 19.0, 22.0, 23.0, 24.0, 26.0);
	cout << "mat33:" << endl << mat33 << endl;
	cout << "det(mat33) = " << mat33.det() << "  (right result = -657)" << endl; 
	cout << "inv(mat33) = " << endl << mat33.inv(mat33.det()) << endl;
	cout << "mat44:" << endl << mat44 << endl;
	cout << "det(mat44) = " << mat44.det() << "  (right result = -320)" << endl; 
	cout << "inv(mat4) = " << endl << mat44.inv(mat44.det()) << endl;
	cout << "mat55:" << endl << mat55 << endl;
	cout << "det(mat55) = " << mat55.det() << "  (right result = 1.001E4)" << endl; 
	cout << "inv(mat55) = " << endl << mat55.inv(mat55.det()) << endl;

	cout << "Checking of the Mat_Band_*m class:" << endl;
	Mat_Band_cm M_U(5,3);
	M_U[1][1] = 10;
	M_U[1][2] = 9;
	M_U[1][3] = 8;
	M_U[2][2] = 20;
	M_U[2][3] = 5;
	M_U[2][4] = 4;
	//M_U[3][2] = 5;
	M_U[3][3] = 30;
	M_U[3][4] = 6;
	M_U[3][5] = 3;
	M_U[4][4] = 40;
	M_U[4][5] = 1;
	M_U[5][5] = 50;

	cout << "upper defined M_U:" << endl << M_U << endl;
	Mat_Band_cm M_L(5,3);
	M_U[1][1] = 10;
	M_U[2][1] = 9;
	M_U[3][1] = 8;
	M_U[2][2] = 20;
	M_U[3][2] = 5;
	M_U[4][2] = 4;
	M_U[3][3] = 30;
	M_U[4][3] = 6;
	M_U[5][3] = 3;
	M_U[4][4] = 40;
	M_U[5][4] = 1;
	M_U[5][5] = 50;
	cout << "lower defined M_L:" << endl << M_U << endl;

	Mat<3,6> MB(1.0f,2.0f,3.0f,4.0f,5.0f,6.0f,7.0f,8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f,16.0f,17.0f,18.0f);
	Mat<3,3> D(11.0f,12.0f,13.0f,12.0f,22.0f,23.0f,13.0f,23.0f,33.0f);
	cout << "MB:" << endl << MB;
	cout << "D:" << endl << D;
	cout << "MB.t*D*MB:" << endl << MB.transpose() * D* MB << endl;
	Mat<6,3> tmp = MB.transpose()*D;
	cout << "tmp = MB.t*D;tmp*MB" << endl << tmp*MB << endl;
	Vec<3> ve(1.0f,2.0f,3.0f);
	ve+=Vec<3>(10.0f, 20.0f, 30.0f);
	cout << "ve:" << ve << endl;

	dMat dmat1(3,4);
	cout << "dmat1:" << endl << dmat1 << endl; //выдает фигню. надо инициализировать массив нулями
	dMat dmat2(3,3,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0);
	cout << "dmat2:" << endl << dmat2 << endl;
	cout << "mat33:" << endl << mat33 << endl;;
	dmat1.cpMat<3,3>(mat33);
	cout << "dmat1.cpMat<3,3>(mat33):" << endl << dmat1 << endl;
	mat33 = dmat2.toMat<3,3>();
	cout << "mat33 = dmat2.toMat<3,3>()" << endl << mat33 << endl;
	mat33 = test_func_dmat().toMat<3,3>();
	cout << "mat33 = test_func_dmat().toMat<3,3>()" << endl << mat33 << endl;
	cout << "dmat1:" << endl << dmat1 << endl;
	cout << "dmat2:" << endl << dmat2 << endl;
	dmat1 = dmat2;
	cout << "dmat1 = dmat2:" << endl << dmat1 << endl;
}

dMat test_func_dmat()
{
	return dMat(3,3,1.1,4.4,7.7,5.5,2.2,8.8,6.6,9.9,3.3);
}

void MKL_test () {
	int n = 5;
	int kd = 2;
	int ldab = kd+1;
	int info;
	int nrhs = 1;
	int ldb = n;
	int *ipv = new int[n];
	Mat_Band_cm M_U(5,3);
	M_U[1][1] = 10;
	M_U[1][2] = 9;
	M_U[1][3] = 8;
	M_U[2][2] = 20;
	M_U[2][3] = 5;
	M_U[2][4] = 4;
	M_U[3][3] = 30;
	M_U[3][4] = 6;
	M_U[3][5] = 3;
	M_U[4][4] = 40;
	M_U[4][5] = 1;
	M_U[5][5] = 50;
	Vec_Long mat(15);
	mat[1]=10;
	mat[2]=9;
	mat[3]=20;
	mat[4]=8;
	mat[5]=5;
	mat[6]=30;
	mat[7]=0;
	mat[8]=4;
	mat[9]=6;
	mat[10]=40;
	mat[11]=0;
	mat[12]=0;
	mat[13]=3;
	mat[14]=1;
	mat[15]=50;

	Vec_Long B(5);
	B[1] = 1;
	B[2] = 0;
	B[3] = 10;
	B[4] = 15;
	B[5] = 4;
	//dpbtrf("U",&n,&kd,M_U.Ptr(),&ldab,&info);
	//dpbtrs("U",&n,&kd,&nrhs,M_U.Ptr(),&ldab,B.Ptr(),&ldb,&info);
	//cout << B << endl;
	//результат должен быть
	//9.55252e-005    -0.139571       0.281898        0.345268        0.0561808

	//второй метод
	dsptrf("U", &n, mat.Ptr(), ipv, &info);
	dsptrs("U", &n, &nrhs,  mat.Ptr(), ipv, B.Ptr(), &ldb, &info);
	cout << B << endl;
	delete[] ipv;
}

void Sparse_Mat_test()
{
	int n = 8;
	Sparse_Matrix_SUR spmat(8);
	spmat.start_training();
	//// 7.0  0.0  1.0  0.0  0.0  2.0  7.0  0.0
	//// 0.0 -4.0  8.0  0.0  2.0  0.0  0.0  0.0
	//// 1.0  8.0  1.0  0.0  0.0  0.0  0.0  5.0
	//// 0.0  0.0  0.0  7.0  0.0  0.0  9.0  0.0
	//// 0.0  2.0  0.0  0.0  5.0  1.0  5.0  0.0
	//// 2.0  0.0  0.0  0.0  1.0 -1.0  0.0  5.0
	//// 7.0  0.0  0.0  9.0  5.0  0.0 11.0  0.0
	//// 0.0  0.0  5.0  0.0  0.0  5.0  0.0  5.0
	spmat.add_value(1,1,7.0);
	spmat.add_value(1,3,1.0);
	spmat.add_value(1,6,2.0);
	spmat.add_value(1,7,7.0);
	spmat.add_value(2,5,2.0);
	spmat.add_value(2,3,8.0);
	spmat.add_value(2,2,-4.0);
	spmat.add_value(3,3,1.0);
	spmat.add_value(3,8,5.0);
	spmat.add_value(4,7,9.0);
	spmat.add_value(4,4,7.0);
	spmat.add_value(5,5,5.0);
	spmat.add_value(5,6,1.0);
	spmat.add_value(5,7,5.0);
	spmat.add_value(6,8,5.0);
	spmat.add_value(6,6,-1.0);
	spmat.add_value(7,7,11.0);
	spmat.add_value(8,8,5.0);
	spmat.stop_training();
	//spmat.print_data(); // works correctly
	double b[] = {1.0, 2.0, 4.0, 6.0, 5.0, 3.0, 7.0, 3.0};
	cout << "SUR*vec: " << endl;
	for (int32 i=0; i < 8; i++)
		cout << spmat.mult_vec_i(b,i+1) << endl;

	Sparse_Matrix_R spmat2(8,8);
	spmat2.start_training();
	//// 7.0  0.0  1.0  0.0  0.0  2.0  7.0  0.0
	//// 0.0 -4.0  8.0  0.0  2.0  0.0  0.0  0.0
	//// 1.0  8.0  1.0  0.0  0.0  0.0  0.0  5.0
	//// 0.0  0.0  0.0  7.0  0.0  0.0  9.0  0.0
	//// 0.0  2.0  0.0  0.0  5.0  1.0  5.0  0.0
	//// 2.0  0.0  0.0  0.0  1.0 -1.0  0.0  5.0
	//// 7.0  0.0  0.0  9.0  5.0  0.0 11.0  0.0
	//// 0.0  0.0  5.0  0.0  0.0  5.0  0.0  5.0
	spmat2.add_value(1,1,7.0);
	spmat2.add_value(3,1,1.0);
	spmat2.add_value(1,6,2.0);
	spmat2.add_value(1,7,7.0);
	spmat2.add_value(2,5,2.0);
	spmat2.add_value(2,3,8.0);
	spmat2.add_value(2,2,-4.0);
	spmat2.add_value(3,3,1.0);
	spmat2.add_value(3,8,5.0);
	spmat2.add_value(4,7,9.0);
	spmat2.add_value(4,4,7.0);
	spmat2.add_value(5,5,5.0);
	spmat2.add_value(5,6,1.0);
	spmat2.add_value(5,7,5.0);
	spmat2.add_value(6,8,5.0);
	spmat2.add_value(6,6,-1.0);
	spmat2.add_value(7,7,11.0);
	spmat2.add_value(8,8,5.0);
	spmat2.add_value(7,1,7.0);
	spmat2.stop_training();
	cout << "R*vec: " << endl;
	for (int32 i=0; i < 8; i++)
		cout << spmat2.mult_vec_i(b,i+1) << endl;
	cout << "R^T*vec: " << endl;
	for (int32 i=0; i < 8; i++)
		cout << spmat2.transpose_mult_vec_i(b,i+1) << endl;
	//double x[8];
	//int nrhs = 1;
	//int mtype = -2; //real symmetric matrix
	//void *pt[64]; /*Internal solver memory pointer pt,*/
	//int iparm[64];
	//int maxfct, mnum, phase, error, msglvl;
	//int i;
	//for (i=0;i<64;i++) iparm[i]=0;
	//iparm[0] = 1; //no solver default
	//iparm[1] = 2; //fill-in reordering from meris
	//iparm[2] = MKL_Get_Max_Threads();
	//iparm[3] = 0; //no iterative-direct algorithm
	//iparm[4] = 0; //no user fill-in reducing permutation
	//iparm[5] = 0; //write solution into x
	//iparm[6] = 16; //default logical fortran unit number for output
	//iparm[7] = 2; //max numbers of iterative refinement steps
	//iparm[9] = 13; //pertrub the pivor elements with 1e-13
	//iparm[10] = 1; //use nonsymmetric permutation  and scaling MPS
	//iparm[13]=0; //output: number of perturbed pivots
	//iparm[17]=-1; //output: number of nonzeros in the factor LU
	//iparm[18]=-1; //output: MFLOPS for LU factorization
	//iparm[19] = 0; //output: number of CG Iterations

	//maxfct = 1; //maximum number of numerical factorizations
	//mnum = 1; //wich factorization to use
	//msglvl = 0; //dont print statistical information in file
	//error = 0; //init error flag
	//for (i=0;i<64;i++) pt[i]=0;

	////Reordering and Symbolic Factorization. This step also allocate 
	//// all memory that is necessary for the factorization
	//phase = 11;
	//PARDISO(pt, &maxfct, &mnum, &mtype,&phase,
	//	&n, spmat.get_values_array(), (int*) spmat.get_iofeir_array(), (int*) spmat.get_columns_array(), NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
	//if (error != 0)
	//{
	//	printf("\nERROR during symbolic factorization: %d", error);
	//	exit(1);
	//}

	//printf("\nReordering completed...");
	//printf("\nNumber of nonzeros in factors=%d", iparm[17]);
	//printf("\nNumber of factorization MFLOPS=%d", iparm[18]);

	//// Numerical factorization

	//phase = 22;
	//PARDISO(pt, &maxfct, &mnum, &mtype,&phase,
	//		&n, spmat.get_values_array(), (int*) spmat.get_iofeir_array(), (int*) spmat.get_columns_array(), NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
	//if (error != 0)
	//{
	//	printf("\nERROR during numerical factorization: %d", error);
	//	exit(2);
	//}
	//printf("\nFactorization completed...");

	//////Back substitution and iterative refiment

	//phase = 33;
	//iparm[7] = 2; //max num of iterative refiment
	//////set right hand b[i]
	//PARDISO(pt, &maxfct, &mnum, &mtype,&phase,
	//		&n, spmat.get_values_array(), (int*) spmat.get_iofeir_array(), (int*) spmat.get_columns_array(), NULL, &nrhs, iparm, &msglvl, b, x, &error);
	//if (error != 0)
	//{
	//	printf("\nERROR during solution: %d", error);
	//	exit(3);
	//}
	//printf("\nSolve completed...");
	//printf("\nThe solution of the system is:\n");
	////output
	//for (i=0;i<n;i++) cout << "  " << x[i] << endl;

	//////Termination and release of memory
	//phase = 1;

	//PARDISO(pt, &maxfct, &mnum, &mtype,&phase,
	//		&n, NULL, (int*) spmat.get_iofeir_array(), (int*) spmat.get_columns_array(), NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
}

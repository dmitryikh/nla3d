#include "sys.h"
#include <mkl.h>
#include "math\SparseMatrix.h"

using namespace std;
using namespace nla3d::math;
//TODO: undone! (see commented lines below)
int main()
{
	int n = 8;
	SparseSymmetricMatrix spmat(8);
	spmat.startTraining();
	//// 7.0  0.0  1.0  0.0  0.0  2.0  7.0  0.0
	//// 0.0 -4.0  8.0  0.0  2.0  0.0  0.0  0.0
	//// 1.0  8.0  1.0  0.0  0.0  0.0  0.0  5.0
	//// 0.0  0.0  0.0  7.0  0.0  0.0  9.0  0.0
	//// 0.0  2.0  0.0  0.0  5.0  1.0  5.0  0.0
	//// 2.0  0.0  0.0  0.0  1.0 -1.0  0.0  5.0
	//// 7.0  0.0  0.0  9.0  5.0  0.0 11.0  0.0
	//// 0.0  0.0  5.0  0.0  0.0  5.0  0.0  5.0
	spmat.addValue(1,1,7.0);
	spmat.addValue(1,3,1.0);
	spmat.addValue(1,6,2.0);
	spmat.addValue(1,7,7.0);
	spmat.addValue(2,5,2.0);
	spmat.addValue(2,3,8.0);
	spmat.addValue(2,2,-4.0);
	spmat.addValue(3,3,1.0);
	spmat.addValue(3,8,5.0);
	spmat.addValue(4,7,9.0);
	spmat.addValue(4,4,7.0);
	spmat.addValue(5,5,5.0);
	spmat.addValue(5,6,1.0);
	spmat.addValue(5,7,5.0);
	spmat.addValue(6,8,5.0);
	spmat.addValue(6,6,-1.0);
	spmat.addValue(7,7,11.0);
	spmat.addValue(8,8,5.0);
	spmat.stopTraining();
	//spmat.print_data(); // works correctly
	double b[] = {1.0, 2.0, 4.0, 6.0, 5.0, 3.0, 7.0, 3.0};
	cout << "SUR*vec: " << endl;
	for (int32 i=0; i < 8; i++)
		cout << spmat.mult_vec_i(b,i+1) << endl;

	SparseMatrix spmat2(8,8);
	spmat2.startTraining();
	//// 7.0  0.0  1.0  0.0  0.0  2.0  7.0  0.0
	//// 0.0 -4.0  8.0  0.0  2.0  0.0  0.0  0.0
	//// 1.0  8.0  1.0  0.0  0.0  0.0  0.0  5.0
	//// 0.0  0.0  0.0  7.0  0.0  0.0  9.0  0.0
	//// 0.0  2.0  0.0  0.0  5.0  1.0  5.0  0.0
	//// 2.0  0.0  0.0  0.0  1.0 -1.0  0.0  5.0
	//// 7.0  0.0  0.0  9.0  5.0  0.0 11.0  0.0
	//// 0.0  0.0  5.0  0.0  0.0  5.0  0.0  5.0
	spmat2.addValue(1,1,7.0);
	spmat2.addValue(3,1,1.0);
	spmat2.addValue(1,6,2.0);
	spmat2.addValue(1,7,7.0);
	spmat2.addValue(2,5,2.0);
	spmat2.addValue(2,3,8.0);
	spmat2.addValue(2,2,-4.0);
	spmat2.addValue(3,3,1.0);
	spmat2.addValue(3,8,5.0);
	spmat2.addValue(4,7,9.0);
	spmat2.addValue(4,4,7.0);
	spmat2.addValue(5,5,5.0);
	spmat2.addValue(5,6,1.0);
	spmat2.addValue(5,7,5.0);
	spmat2.addValue(6,8,5.0);
	spmat2.addValue(6,6,-1.0);
	spmat2.addValue(7,7,11.0);
	spmat2.addValue(8,8,5.0);
	spmat2.addValue(7,1,7.0);
	spmat2.stopTraining();
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

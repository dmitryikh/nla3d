#include "sys.h"
#include "math/SparseMatrix.h"

using namespace std;
using namespace nla3d::math;


int main() {
  cout << "SparceMatrix(3,3) with no elements" << endl;
  {
    SparseMatrix spmat(3, 3);
    spmat.startTraining();
    // add nothing
    spmat.stopTraining();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparceMatrix(3, 0) with no elements" << endl;
  {
    SparseMatrix spmat(3, 0);
    spmat.startTraining();
    // add nothing
    spmat.stopTraining();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparceMatrix(0, 3) with no elements" << endl;
  {
    SparseMatrix spmat(0, 3);
    spmat.startTraining();
    // add nothing
    spmat.stopTraining();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparseSymmetricMatrix(3,3) with no elements" << endl;
  {
    SparseSymmetricMatrix spmat(3);
    spmat.startTraining();
    // add nothing
    spmat.stopTraining();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparseSymmetricMatrix(3, 0) with no elements" << endl;
  {
    SparseSymmetricMatrix spmat(0);
    spmat.startTraining();
    // add nothing
    spmat.stopTraining();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }





	// int n = 8;
	// SparseSymmetricMatrix spmat(8);
	// spmat.startTraining();
	// //// 7.0  0.0  1.0  0.0  0.0  2.0  7.0  0.0
	// //// 0.0 -4.0  8.0  0.0  2.0  0.0  0.0  0.0
	// //// 1.0  8.0  1.0  0.0  0.0  0.0  0.0  5.0
	// //// 0.0  0.0  0.0  7.0  0.0  0.0  9.0  0.0
	// //// 0.0  2.0  0.0  0.0  5.0  1.0  5.0  0.0
	// //// 2.0  0.0  0.0  0.0  1.0 -1.0  0.0  5.0
	// //// 7.0  0.0  0.0  9.0  5.0  0.0 11.0  0.0
	// //// 0.0  0.0  5.0  0.0  0.0  5.0  0.0  5.0
	// spmat.addValue(1,1,7.0);
	// spmat.addValue(1,3,1.0);
	// spmat.addValue(1,6,2.0);
	// spmat.addValue(1,7,7.0);
	// spmat.addValue(2,5,2.0);
	// spmat.addValue(2,3,8.0);
	// spmat.addValue(2,2,-4.0);
	// spmat.addValue(3,3,1.0);
	// spmat.addValue(3,8,5.0);
	// spmat.addValue(4,7,9.0);
	// spmat.addValue(4,4,7.0);
	// spmat.addValue(5,5,5.0);
	// spmat.addValue(5,6,1.0);
	// spmat.addValue(5,7,5.0);
	// spmat.addValue(6,8,5.0);
	// spmat.addValue(6,6,-1.0);
	// spmat.addValue(7,7,11.0);
	// spmat.addValue(8,8,5.0);
	// spmat.stopTraining();
	// //spmat.print_data(); // works correctly
	// double b[] = {1.0, 2.0, 4.0, 6.0, 5.0, 3.0, 7.0, 3.0};
	// cout << "SUR*vec: " << endl;
	// for (int32 i=0; i < 8; i++)
	// 	cout << spmat.mult_vec_i(b,i+1) << endl;

	// SparseMatrix spmat2(8,8);
	// spmat2.startTraining();
	// //// 7.0  0.0  1.0  0.0  0.0  2.0  7.0  0.0
	// //// 0.0 -4.0  8.0  0.0  2.0  0.0  0.0  0.0
	// //// 1.0  8.0  1.0  0.0  0.0  0.0  0.0  5.0
	// //// 0.0  0.0  0.0  7.0  0.0  0.0  9.0  0.0
	// //// 0.0  2.0  0.0  0.0  5.0  1.0  5.0  0.0
	// //// 2.0  0.0  0.0  0.0  1.0 -1.0  0.0  5.0
	// //// 7.0  0.0  0.0  9.0  5.0  0.0 11.0  0.0
	// //// 0.0  0.0  5.0  0.0  0.0  5.0  0.0  5.0
	// spmat2.addValue(1,1,7.0);
	// spmat2.addValue(3,1,1.0);
	// spmat2.addValue(1,6,2.0);
	// spmat2.addValue(1,7,7.0);
	// spmat2.addValue(2,5,2.0);
	// spmat2.addValue(2,3,8.0);
	// spmat2.addValue(2,2,-4.0);
	// spmat2.addValue(3,3,1.0);
	// spmat2.addValue(3,8,5.0);
	// spmat2.addValue(4,7,9.0);
	// spmat2.addValue(4,4,7.0);
	// spmat2.addValue(5,5,5.0);
	// spmat2.addValue(5,6,1.0);
	// spmat2.addValue(5,7,5.0);
	// spmat2.addValue(6,8,5.0);
	// spmat2.addValue(6,6,-1.0);
	// spmat2.addValue(7,7,11.0);
	// spmat2.addValue(8,8,5.0);
	// spmat2.addValue(7,1,7.0);
	// spmat2.stopTraining();
	// cout << "R*vec: " << endl;
	// for (int32 i=0; i < 8; i++)
	// 	cout << spmat2.mult_vec_i(b,i+1) << endl;
	// cout << "R^T*vec: " << endl;
	// for (int32 i=0; i < 8; i++)
	// 	cout << spmat2.transpose_mult_vec_i(b,i+1) << endl;
}

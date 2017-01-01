#include "sys.h"
#include "math/Vec.h"
#include "math/SparseMatrix.h"

using namespace std;
using namespace nla3d::math;


int main() {
  cout << "SparceMatrix(3,3) with no elements" << endl;
  {
    SparseMatrix spmat(3, 3);
    // add nothing
    spmat.compress();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparceMatrix(3, 0) with no elements" << endl;
  {
    SparseMatrix spmat(3, 0);
    // add nothing
    spmat.compress();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparceMatrix(0, 3) with no elements" << endl;
  {
    SparseMatrix spmat(0, 3);
    // add nothing
    spmat.compress();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparseSymmetricMatrix(3,3) with no elements" << endl;
  {
    SparseSymMatrix spmat(3);
    // add nothing
    spmat.compress();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparseSymmetricMatrix(3, 0) with no elements" << endl;
  {
    SparseSymMatrix spmat(0);
    // add nothing
    spmat.compress();

    spmat.printInternalData(cout);
    spmat.print(cout);
  }

  cout << "SparseMatrix(3, 4) with elements" << endl;
  cout << "[ 0 1 0 4]" << endl;
  cout << "[ 0 0 0 0]" << endl;
  cout << "[ 2 8 5 9]" << endl;
  {
    SparseMatrix spmat(3, 4);
    spmat.add(1, 2);
    spmat.add(1, 4);
    spmat.add(3, 1);
    spmat.add(3, 2);
    spmat.add(3, 3);
    spmat.add(3, 4);
    spmat.compress();
    spmat.addValue(1, 2, 1.0);
    spmat.addValue(1, 4, 4.0);
    spmat.addValue(3, 1, 2.0);
    spmat.addValue(3, 2, 8.0);
    spmat.addValue(3, 3, 5.0);
    spmat.addValue(3, 4, 9.0);

    spmat.printInternalData(cout);
    spmat.print(cout);
    // values = {1	4	2	8	5	9	}
    // columns = {2	4	1	2	3	4	}
    // iofeir = {0	2	2	6	}
    CHECK_EQ(spmat.getValuesArray()[0], 1.0);
    CHECK_EQ(spmat.getValuesArray()[1], 4.0);
    CHECK_EQ(spmat.getValuesArray()[2], 2.0);
    CHECK_EQ(spmat.getValuesArray()[3], 8.0);
    CHECK_EQ(spmat.getValuesArray()[4], 5.0);
    CHECK_EQ(spmat.getValuesArray()[5], 9.0);

    CHECK_EQ(spmat.getColumnsArray()[0], 2);
    CHECK_EQ(spmat.getColumnsArray()[1], 4);
    CHECK_EQ(spmat.getColumnsArray()[2], 1);
    CHECK_EQ(spmat.getColumnsArray()[3], 2);
    CHECK_EQ(spmat.getColumnsArray()[4], 3);
    CHECK_EQ(spmat.getColumnsArray()[5], 4);

    CHECK_EQ(spmat.getIofeirArray()[0], 1);
    CHECK_EQ(spmat.getIofeirArray()[1], 3);
    CHECK_EQ(spmat.getIofeirArray()[2], 3);
    CHECK_EQ(spmat.getIofeirArray()[3], 7);
  }


  cout << "SparseSymMatrix(4, 4) with elements" << endl;
  cout << "[ 2 3 0 4]" << endl;
  cout << "[ s 0 8 0]" << endl;
  cout << "[ s s 1 0]" << endl;
  cout << "[ s s s 9]" << endl;
  {
    SparseSymMatrix spmat(4);
    spmat.add(1, 1);
    spmat.add(1, 2);
    spmat.add(1, 4);
    spmat.add(3, 2);
    spmat.add(3, 3);
    spmat.add(4, 4);
    spmat.compress();
    spmat(1, 1) = 2.0;
    spmat(2, 1) = 3.0;
    spmat(4, 1) = 4.0;
    spmat(2, 3) = 8.0;
    spmat(3, 3) = 1.0;
    spmat(4, 4) = 9.0;

    spmat.printInternalData(cout);
    spmat.print(cout);
    // values = {2	3	4	8	1	9	}
    // columns = {1	2	4	3	3	4	}
    // iofeir = {0	3	4	5	6	}
    CHECK_EQ(spmat.getValuesArray()[0], 2.0);
    CHECK_EQ(spmat.getValuesArray()[1], 3.0);
    CHECK_EQ(spmat.getValuesArray()[2], 4.0);
    CHECK_EQ(spmat.getValuesArray()[3], 0.0);
    CHECK_EQ(spmat.getValuesArray()[4], 8.0);
    CHECK_EQ(spmat.getValuesArray()[5], 1.0);
    CHECK_EQ(spmat.getValuesArray()[6], 9.0);

    CHECK_EQ(spmat.getColumnsArray()[0], 1);
    CHECK_EQ(spmat.getColumnsArray()[1], 2);
    CHECK_EQ(spmat.getColumnsArray()[2], 4);
    CHECK_EQ(spmat.getColumnsArray()[3], 2);
    CHECK_EQ(spmat.getColumnsArray()[4], 3);
    CHECK_EQ(spmat.getColumnsArray()[5], 3);
    CHECK_EQ(spmat.getColumnsArray()[6], 4);

    CHECK_EQ(spmat.getIofeirArray()[0], 1);
    CHECK_EQ(spmat.getIofeirArray()[1], 4);
    CHECK_EQ(spmat.getIofeirArray()[2], 6);
    CHECK_EQ(spmat.getIofeirArray()[3], 7);
    CHECK_EQ(spmat.getIofeirArray()[4], 8);
  }

  cout << "Create two matrices with the same SparsityInfo" << endl;
  cout << "A = " << endl;
  cout << "[ 0 0 1 0]" << endl;
  cout << "[ 1 0 1 0]" << endl;
  cout << "[ 0 0 0 0]" << endl;
  cout << "[ 0 0 1 0]" << endl;

  cout << "B = " << endl;
  cout << "[ 2 0 0 2]" << endl;
  cout << "[ 0 0 2 0]" << endl;
  cout << "[ 0 2 0 0]" << endl;
  cout << "[ 0 0 0 0]" << endl;
  {
    SparseMatrix A(4,4);
    A.add(1, 3);
    A.add(2, 1);
    {
      // share SparsityInfo
      SparseMatrix B(A.getSparsityInfo());
      
      A.add(2, 3);
      A.add(4, 3);

      B.add(1, 1);
      B.add(1, 4);
      B.add(2, 3);
      B.add(3, 2);

      A.compress();
      B.compress();

      A(1, 3) = 1.0;
      A(2, 1) = 1.0;
      A(2, 3) = 1.0;
      A(4, 3) = 1.0;

      B(1, 1) = 2.0;
      B(1, 4) = 2.0;
      B(2, 3) = 2.0;
      B(3, 2) = 2.0;

      cout << "B structure:" << endl;
      B.printInternalData(cout);
      B.print(cout);
      // values = {2	0	2	0	2	2	0	}
      // columns = {1	3	4	1	3	2	3	}
      // iofeir = {0	3	5	6	7	}
      CHECK_EQ(B.getValuesArray()[0], 2.0);
      CHECK_EQ(B.getValuesArray()[1], 0.0);
      CHECK_EQ(B.getValuesArray()[2], 2.0);
      CHECK_EQ(B.getValuesArray()[3], 0.0);
      CHECK_EQ(B.getValuesArray()[4], 2.0);
      CHECK_EQ(B.getValuesArray()[5], 2.0);
      CHECK_EQ(B.getValuesArray()[6], 0.0);

      CHECK_EQ(B.getColumnsArray()[0], 1);
      CHECK_EQ(B.getColumnsArray()[1], 3);
      CHECK_EQ(B.getColumnsArray()[2], 4);
      CHECK_EQ(B.getColumnsArray()[3], 1);
      CHECK_EQ(B.getColumnsArray()[4], 3);
      CHECK_EQ(B.getColumnsArray()[5], 2);
      CHECK_EQ(B.getColumnsArray()[6], 3);

      CHECK_EQ(B.getIofeirArray()[0], 1);
      CHECK_EQ(B.getIofeirArray()[1], 4);
      CHECK_EQ(B.getIofeirArray()[2], 6);
      CHECK_EQ(B.getIofeirArray()[3], 7);
      CHECK_EQ(B.getIofeirArray()[4], 8);
      // destroy B here
    }

    cout << "A structure:" << endl;
    A.printInternalData(cout);
    A.print(cout);
    // values = {0	1	0	1	1	0	1	}
    // columns = {1	3	4	1	3	2	3	}
    // iofeir = {0	3	5	6	7	}
    CHECK_EQ(A.getValuesArray()[0], 0.0);
    CHECK_EQ(A.getValuesArray()[1], 1.0);
    CHECK_EQ(A.getValuesArray()[2], 0.0);
    CHECK_EQ(A.getValuesArray()[3], 1.0);
    CHECK_EQ(A.getValuesArray()[4], 1.0);
    CHECK_EQ(A.getValuesArray()[5], 0.0);
    CHECK_EQ(A.getValuesArray()[6], 1.0);

    CHECK_EQ(A.getColumnsArray()[0], 1);
    CHECK_EQ(A.getColumnsArray()[1], 3);
    CHECK_EQ(A.getColumnsArray()[2], 4);
    CHECK_EQ(A.getColumnsArray()[3], 1);
    CHECK_EQ(A.getColumnsArray()[4], 3);
    CHECK_EQ(A.getColumnsArray()[5], 2);
    CHECK_EQ(A.getColumnsArray()[6], 3);

    CHECK_EQ(A.getIofeirArray()[0], 1);
    CHECK_EQ(A.getIofeirArray()[1], 4);
    CHECK_EQ(A.getIofeirArray()[2], 6);
    CHECK_EQ(A.getIofeirArray()[3], 7);
    CHECK_EQ(A.getIofeirArray()[4], 8);

  }

  cout << "SparseMatrix(3, 4) math" << endl;
  cout << "    [ 0 1 0 4]" << endl;
  cout << "A = [ 0 0 0 0]" << endl;
  cout << "    [ 2 8 5 9]" << endl;
  {
    SparseMatrix A(3, 4);
    A.add(1, 2);
    A.add(1, 4);
    A.add(3, 1);
    A.add(3, 2);
    A.add(3, 3);
    A.add(3, 4);
    A.compress();
    A(1, 2) = 1.0;
    A(1, 4) = 4.0;
    A(3, 1) = 2.0;
    A(3, 2) = 8.0;
    A(3, 3) = 5.0;
    A(3, 4) = 9.0;

    Vec<4> b1 = {1.0, 2.0, 3.0, 4.0};
    Vec<3> b2 = {3.0, 2.0, 1.0};
    Vec<4> res;

    cout << "b1 = " << b1 << endl;
    cout << "b2 = " << b2 << endl;

    cout << "A.mult_vec_i(b1, 1) = " << A.mult_vec_i(b1.ptr(), 1) << endl;
    cout << "A.mult_vec_i(b1, 2) = " << A.mult_vec_i(b1.ptr(), 2) << endl;
    cout << "A.mult_vec_i(b1, 3) = " << A.mult_vec_i(b1.ptr(), 3) << endl;

    CHECK_EQ(A.mult_vec_i(b1.ptr(), 1), 18.0);
    CHECK_EQ(A.mult_vec_i(b1.ptr(), 2),  0.0);
    CHECK_EQ(A.mult_vec_i(b1.ptr(), 3), 69.0);

    cout << "A.transpose_mult_vec_i(b2, 1) = " << A.transpose_mult_vec_i(b2.ptr(), 1) << endl;
    cout << "A.transpose_mult_vec_i(b2, 2) = " << A.transpose_mult_vec_i(b2.ptr(), 2) << endl;
    cout << "A.transpose_mult_vec_i(b2, 3) = " << A.transpose_mult_vec_i(b2.ptr(), 3) << endl;
    cout << "A.transpose_mult_vec_i(b2, 4) = " << A.transpose_mult_vec_i(b2.ptr(), 4) << endl;

    CHECK_EQ(A.transpose_mult_vec_i(b2.ptr(), 1),  2.0);
    CHECK_EQ(A.transpose_mult_vec_i(b2.ptr(), 2), 11.0);
    CHECK_EQ(A.transpose_mult_vec_i(b2.ptr(), 3),  5.0);
    CHECK_EQ(A.transpose_mult_vec_i(b2.ptr(), 4), 21.0);

    A.transpose_mult_vec(b2.ptr(), res.ptr());
    cout << "A.transpose_mult_vec(b2, res); res = " << res << endl;

    CHECK_EQ(res[0],  2.0);
    CHECK_EQ(res[1], 11.0);
    CHECK_EQ(res[2],  5.0);
    CHECK_EQ(res[3], 21.0);
  }


  cout << "SparseSymMatrix(8, 8) math" << endl;
  cout << "    [7.0  0.0  1.0  0.0  0.0  2.0  7.0  0.0]" << endl;
  cout << "    [ s  -4.0  8.0  0.0  2.0  0.0  0.0  0.0]" << endl;
  cout << "    [ s    s   1.0  0.0  0.0  0.0  0.0  5.0]" << endl;
  cout << "A = [ s    s    s   7.0  0.0  0.0  9.0  0.0]" << endl;
  cout << "    [ s    s    s    s   5.0  1.0  5.0  0.0]" << endl;
  cout << "    [ s    s    s    s    s  -1.0  0.0  5.0]" << endl;
  cout << "    [ s    s    s    s    s    s  11.0  0.0]" << endl;
  cout << "    [ s    s    s    s    s    s    s   5.0]" << endl;
  // Expression for Python to check math by numpy
	// A=np.array([[7.0, 0.0, 1.0, 0.0, 0.0, 2.0, 7.0, 0.0],
	//             [0.0,-4.0, 8.0, 0.0, 2.0, 0.0, 0.0, 0.0],
	//             [1.0, 8.0, 1.0, 0.0, 0.0, 0.0, 0.0, 5.0],
	//             [0.0, 0.0, 0.0, 7.0, 0.0, 0.0, 9.0, 0.0],
	//             [0.0, 2.0, 0.0, 0.0, 5.0, 1.0, 5.0, 0.0],
	//             [2.0, 0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 5.0],
	//             [7.0, 0.0, 0.0, 9.0, 5.0, 0.0,11.0, 0.0],
	//             [0.0, 0.0, 5.0, 0.0, 0.0, 5.0, 0.0, 5.0]])
  {
    SparseSymMatrix A(8);

	  A.add(1, 1);
	  A.add(1, 3);
	  A.add(1, 6);
	  A.add(1, 7);
	  A.add(2, 5);
	  A.add(2, 3);
	  A.add(2, 2);
	  A.add(3, 3);
	  A.add(3, 8);
	  A.add(4, 7);
	  A.add(4, 4);
	  A.add(5, 5);
	  A.add(5, 6);
	  A.add(5, 7);
	  A.add(6, 8);
	  A.add(6, 6);
	  A.add(7, 7);
	  A.add(8, 8);
    
    A.compress();

	  A(1, 1) =  7.0;
	  A(1, 3) =  1.0;
	  A(1, 6) =  2.0;
	  A(1, 7) =  7.0;
	  A(2, 5) =  2.0;
	  A(2, 3) =  8.0;
	  A(2, 2) = -4.0;
	  A(3, 3) =  1.0;
	  A(3, 8) =  5.0;
	  A(4, 7) =  9.0;
	  A(4, 4) =  7.0;
	  A(5, 5) =  5.0;
	  A(5, 6) =  1.0;
	  A(5, 7) =  5.0;
	  A(6, 8) =  5.0;
	  A(6, 6) = -1.0;
	  A(7, 7) = 11.0;
	  A(8, 8) =  5.0;

	  Vec<8> b = {1.0, 2.0, 4.0, 6.0, 5.0, 3.0, 7.0, 3.0};

    cout << "b = " << b << endl;

    cout << "A.mult_vec_i(b, 1) = " << A.mult_vec_i(b.ptr(), 1) << endl;
    cout << "A.mult_vec_i(b, 2) = " << A.mult_vec_i(b.ptr(), 2) << endl;
    cout << "A.mult_vec_i(b, 3) = " << A.mult_vec_i(b.ptr(), 3) << endl;
    cout << "A.mult_vec_i(b, 4) = " << A.mult_vec_i(b.ptr(), 4) << endl;
    cout << "A.mult_vec_i(b, 5) = " << A.mult_vec_i(b.ptr(), 5) << endl;
    cout << "A.mult_vec_i(b, 6) = " << A.mult_vec_i(b.ptr(), 6) << endl;
    cout << "A.mult_vec_i(b, 7) = " << A.mult_vec_i(b.ptr(), 7) << endl;
    cout << "A.mult_vec_i(b, 8) = " << A.mult_vec_i(b.ptr(), 8) << endl;

    CHECK_EQ(A.mult_vec_i(b.ptr(), 1), 66.0);
    CHECK_EQ(A.mult_vec_i(b.ptr(), 2), 34.0);
    CHECK_EQ(A.mult_vec_i(b.ptr(), 3), 36.0);
    CHECK_EQ(A.mult_vec_i(b.ptr(), 4),105.0);
    CHECK_EQ(A.mult_vec_i(b.ptr(), 5), 67.0);
    CHECK_EQ(A.mult_vec_i(b.ptr(), 6), 19.0);
    CHECK_EQ(A.mult_vec_i(b.ptr(), 7),163.0);
    CHECK_EQ(A.mult_vec_i(b.ptr(), 8), 50.0);
  }

}

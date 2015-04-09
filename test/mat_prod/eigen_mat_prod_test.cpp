#include <fstream>
#include "sys.h"
#include "math/Mat.h"
#include <Eigen/Dense>

using namespace std;
using namespace nla3d;
using namespace nla3d::math;

const double eps = 0.000001;
uint16 tn = 100;
char* dir;

template<typename M1, typename M2>
bool eigen_compare(const Eigen::MatrixBase<M1>& ref, const Eigen::MatrixBase<M2>& res, double eps = 0.001) {
  assert(ref.rows() == res.rows());
  assert(ref.cols() == res.cols());
  for (uint16 i = 0; i < ref.rows(); i++) {
    for (uint16 j = 0; j < ref.cols(); j++) {
      if (fabs(ref(i,j) - res(i,j)) > eps) {
        LOG(ERROR) << "Matrices are different";
        exit(1);
        return false;
      }
    }
  }
  return true;
}

bool test_matBVprod () {
	std::ifstream in;
	const uint16 _N = 24;
	const uint16 _M = 9;
	char filename[100];
	Mat2<_N,_M> B;
	Vec<_M> V;
	Vec<_N> Rf;
	Vec<_N> R;

  Eigen::MatrixXd e_B;
  Eigen::VectorXd e_V;
  Eigen::VectorXd e_R;

  sprintf_s(filename,100,"%s/matBVprod_%02d%02d",dir,_N,_M);
  in.open(filename);
	for (uint16 gg = 1; gg <= tn; gg++) {
		B.zeros();
		V.zeros();
		R.zeros();
		Rf.zeros();

		B.simple_read(in);
		V.simple_read(in);
		Rf.simple_read(in);

		matBVprod(B, V, 1.0, R);

    e_B = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (B.ptr(), _N, _M);
    e_V = Eigen::Map<Eigen::VectorXd> (V.ptr(), _M, 1);

    e_R = e_B*e_V;

		if (!R.compare(Rf, eps)) {
			LOG(ERROR) << "test_matBVprod: n = " << gg;
      exit(1);
		}
    eigen_compare(Eigen::Map<Eigen::MatrixXd> (Rf.ptr(),_N,1),e_R);
		DLOG(DEBUG) << "test_matBVprod: case " << gg << " checked successfuly!";
	}
  in.close();
	return true;
}

bool test_matBTVprod () {
	std::ifstream in;
	const uint16 _N = 24;
	const uint16 _M = 9;
	char filename[100];
	Mat2<_N,_M> B;
	Vec<_N> V;
	Vec<_M> Rf;
	Vec<_M> R;

  Eigen::MatrixXd e_B;
  Eigen::VectorXd e_V;
  Eigen::VectorXd e_R;

  sprintf_s(filename,100, "%s/matBTVprod_%02d%02d", dir, _N, _M);
  in.open(filename);
	for (uint16 gg = 1; gg <= tn; gg++) {
		B.zeros();
		V.zeros();
		R.zeros();
		Rf.zeros();

		B.simple_read(in);
		V.simple_read(in);
		Rf.simple_read(in);

		matBTVprod(B, V, 1.0, R);

    e_B = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (B.ptr(), _N, _M);
    e_V = Eigen::Map<Eigen::VectorXd> (V.ptr(), _N, 1);

    e_R = e_B.transpose()*e_V;

		if (!R.compare(Rf, eps)) {
			LOG(ERROR) << "test_matBTVprod: n = " << gg;
      exit(1);
		}
    eigen_compare(Eigen::Map<Eigen::MatrixXd > (Rf.ptr(),_M,1), e_R);
		DLOG(DEBUG) << "test_matBTVprod: case " << gg << " checked successfuly!";
	}
  in.close();
	return true;
}

bool test_matABprod () {
  std::ifstream in;
	const uint16 _N = 24;
	const uint16 _M = 9;
	const uint16 _M2= 12;
	char filename[100];
	Mat2<_N,_M>  A;
	Mat2<_M,_M2> B;
	Mat2<_N,_M2> R;
	Mat2<_N,_M2> Rf;
  Eigen::MatrixXd e_A, e_B, e_R;
  sprintf_s(filename,100, "%s/matABprod_%02d%02d%02d", dir, _N, _M, _M2);
  in.open(filename);
	for (uint16 gg = 1; gg <= tn; gg++) {
		A.zeros();
		B.zeros();
		R.zeros();
		Rf.zeros();

		A.simple_read(in);
		B.simple_read(in);
		Rf.simple_read(in);

		matABprod(A, B, 1.0, R);

    e_A = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (A.ptr(), _N, _M);
    e_B = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (B.ptr(), _M, _M2);
    e_R = e_A * e_B;
		if (!R.compare(Rf, eps)) {
			LOG(ERROR) << "test_matABprod: n = " << gg;
      exit(1);
		}
    eigen_compare(Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (Rf.ptr(),_N, _M2),
       e_R);
		DLOG(DEBUG) << "test_matABprod: case " << gg << " checked successfuly!";
	}
  in.close();
	return true;
}

bool test_matATBprod () {
  std::ifstream in;
	const uint16 _M = 24;
	const uint16 _N = 9;
	const uint16 _N2= 12;
	char filename[100];
	Mat2<_M,_N>  A;
	Mat2<_M,_N2> B;
	Mat2<_N,_N2> R;
	Mat2<_N,_N2> Rf;

  Eigen::MatrixXd e_A, e_B, e_R;
  sprintf_s(filename,100, "%s/matATBprod_%02d%02d%02d", dir, _M, _N, _N2);
  in.open(filename);
	for (uint16 gg = 1; gg <= tn; gg++) {
		A.zeros();
		B.zeros();
		R.zeros();
		Rf.zeros();

		A.simple_read(in);
		B.simple_read(in);
		Rf.simple_read(in);

		matATBprod(A, B, 1.0, R);


    e_A = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (A.ptr(), _M, _N);
    e_B = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (B.ptr(), _M, _N2);
    e_R = e_A.transpose() * e_B;
		if (!R.compare(Rf, eps)) {
			LOG(ERROR) << "test_matATBprod: n = " << gg;
      exit(1);
		}
    eigen_compare(Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (Rf.ptr(),_N, _N2),
       e_R);
		DLOG(DEBUG) << "test_matATBprod: case " << gg << " checked successfuly!";
	}
  in.close();
	return true;
}

bool test_matBTDBprod () {
  std::ifstream in;
	const uint16 _M = 24;
	const uint16 _N = 9;
	char filename[100];
	Mat2<_M,_N>  B;
	MatSym<_M> D;
	MatSym<_N> R;
	MatSym<_N> Rf;
  Eigen::MatrixXd e_D, e_B, e_R;
  sprintf_s(filename,100, "%s/matBTDBprod_%02d%02d", dir, _M, _N);
  in.open(filename);
	for (uint16 gg = 1; gg <= tn; gg++) {
		D.zeros();
		B.zeros();
		R.zeros();
		Rf.zeros();

		B.simple_read(in);
		D.simple_read(in);
		Rf.simple_read(in);

		matBTDBprod(B, D, 1.0, R);


    e_D.resize(_M,_M);
    uint16 c = 0;
    for (uint16 i = 0; i < _M; i++) {
      for (uint16 j = i; j < _M; j++) {
        e_D(i,j) = D.data[c];
        c++;
      }
    }
    //cout << "e_D = " << endl << e_D << endl;
    //Eigen::MatrixXd tmp = e_D.selfadjointView<Eigen::Upper>();
    //cout << "e_D.symmetric_view = " << endl << tmp;
    e_B = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > (B.ptr(), _M, _N);
    e_R = e_B.transpose() *  e_D.selfadjointView<Eigen::Upper>() * e_B;
    Eigen::MatrixXd e_Rf(_N, _N);
    c = 0;
    for (uint16 i = 0; i < _N; i++) {
      for (uint16 j = i; j < _N; j++) {
        e_Rf(i,j) = Rf.data[c];
        c++;
      }
    }
    e_Rf = e_Rf.selfadjointView<Eigen::Upper>();
    cout << "e_Rf.rows = " << e_Rf.rows() << ",e_Rf.cols = " << e_Rf.cols() << endl;
    cout << "e_R.rows = " << e_R.rows() << ",e_R.cols = " << e_R.cols() << endl;
    eigen_compare(e_Rf, e_R);
		if (!R.compare(Rf, eps)) {
			LOG(ERROR) << "test_matBTDBprod: n = " << gg;
      exit(1);
		}
		DLOG(DEBUG) << "test_matBTDBprod: case " << gg << " checked successfuly!";
	}
  in.close();
	return true;
}

int main (int argc, char* argv[]) {
  char* tmp = getCmdOption(argv, argv + argc, "-dir");
  if (tmp) {
    dir = tmp;
  } else {
    LOG(FATAL) << "You shoud provide directory with test data";
  }

  tmp = getCmdOption(argv, argv + argc, "-num");
  if (tmp) {
    tn = atoi(tmp);
  }
	test_matBVprod();
	test_matBTVprod();
	test_matABprod();
	test_matATBprod();
	test_matBTDBprod();
}

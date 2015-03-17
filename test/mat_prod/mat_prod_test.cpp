#include <fstream>
#include "sys.h"
#include "math/Mat.h"

using namespace std;
using namespace nla3d;
using namespace nla3d::math;

const double eps = 0.000001;
uint16 tn = 100;
char* dir;

bool test_matBVprod ()
{
	std::ifstream in;
	const uint16 _N = 24;
	const uint16 _M = 9;
	char filename[100];
	Mat2<_N,_M> B;
	Vec<_M> V;
	Vec<_N> Rf;
	Vec<_N> R;
  sprintf_s(filename,100,"%s/matBVprod_%02d%02d",dir,_N,_M);
  in.open(filename);
	for (uint16 gg=1;gg<=tn;gg++)
	{
		B.zeros();
		V.zeros();
		R.zeros();
		Rf.zeros();

		B.simple_read(in);
		V.simple_read(in);
		Rf.simple_read(in);

		matBVprod(B, V, 1.0, R);

		if (!R.compare(Rf, eps)) {
			error("test_matBVprod: n=%d", gg);
		}
		debug("test_matBVprod: case %d checked successfuly!", gg);
	}
  in.close();
	return true;
}

bool test_matBTVprod ()
{
	std::ifstream in;
	const uint16 _N = 24;
	const uint16 _M = 9;
	char filename[100];
	Mat2<_N,_M> B;
	Vec<_N> V;
	Vec<_M> Rf;
	Vec<_M> R;
  sprintf_s(filename,100, "%s/matBTVprod_%02d%02d", dir, _N, _M);
  in.open(filename);
	for (uint16 gg=1;gg<=tn;gg++)
	{
		B.zeros();
		V.zeros();
		R.zeros();
		Rf.zeros();

		B.simple_read(in);
		V.simple_read(in);
		Rf.simple_read(in);

		matBTVprod(B, V, 1.0, R);

		if (!R.compare(Rf, eps)) {
			error("test_matBTVprod: n=%d", gg);
		}
		debug("test_matBTVprod: case %d checked successfuly!", gg);
	}
  in.close();
	return true;
}

bool test_matABprod ()
{
std::ifstream in;
	const uint16 _N = 24;
	const uint16 _M = 9;
	const uint16 _M2= 12;
	char filename[100];
	Mat2<_N,_M>  A;
	Mat2<_M,_M2> B;
	Mat2<_N,_M2> R;
	Mat2<_N,_M2> Rf;
  sprintf_s(filename,100, "%s/matABprod_%02d%02d%02d", dir, _N, _M, _M2);
  in.open(filename);
	for (uint16 gg=1;gg<=tn;gg++)
	{
		A.zeros();
		B.zeros();
		R.zeros();
		Rf.zeros();

		A.simple_read(in);
		B.simple_read(in);
		Rf.simple_read(in);

		matABprod(A, B, 1.0, R);

		if (!R.compare(Rf, eps)) {
			error("test_matABprod: n=%d", gg);
		}
		debug("test_matABprod: case %d checked successfuly!", gg);
	}
  in.close();
	return true;
}

bool test_matATBprod ()
{
std::ifstream in;
	const uint16 _M = 24;
	const uint16 _N = 9;
	const uint16 _N2= 12;
	char filename[100];
	Mat2<_M,_N>  A;
	Mat2<_M,_N2> B;
	Mat2<_N,_N2> R;
	Mat2<_N,_N2> Rf;
  sprintf_s(filename,100, "%s/matATBprod_%02d%02d%02d", dir, _M, _N, _N2);
  in.open(filename);
	for (uint16 gg=1;gg<=tn;gg++)
	{
		A.zeros();
		B.zeros();
		R.zeros();
		Rf.zeros();

		A.simple_read(in);
		B.simple_read(in);
		Rf.simple_read(in);

		matATBprod(A, B, 1.0, R);
		if (!R.compare(Rf, eps)) {
			error("test_matATBprod: n=%d", gg);
		}
		debug("test_matATBprod: case %d checked successfuly!", gg);
	}
  in.close();
	return true;
}

bool test_matBTDBprod ()
{
std::ifstream in;
	const uint16 _M = 24;
	const uint16 _N = 9;
	char filename[100];
	Mat2<_M,_N>  B;
	MatSym<_M> D;
	MatSym<_N> R;
	MatSym<_N> Rf;
  sprintf_s(filename,100, "%s/matBTDBprod_%02d%02d", dir, _M, _N);
  in.open(filename);
	for (uint16 gg=1;gg<=tn;gg++)
	{
		D.zeros();
		B.zeros();
		R.zeros();
		Rf.zeros();

		B.simple_read(in);
		D.simple_read(in);
		Rf.simple_read(in);

		matBTDBprod(B, D, 1.0, R);

		if (!R.compare(Rf, eps)) {
			error("test_matBTDBprod: n=%d", gg);
		}
		debug("test_matBTDBprod: case %d checked successfuly!", gg);

	}
  in.close();
	return true;
}

int main (int argc, char* argv[])
{
  char* tmp = getCmdOption(argv, argv + argc, "-dir");
  if (tmp) {
    dir = tmp;
  } else {
    error("You shoud provide directory with test data");
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

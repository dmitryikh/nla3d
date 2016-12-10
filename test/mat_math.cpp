#include "sys.h"
#include "math/Mat.h"

//TODO: not all results are autochecked

using namespace std;
using namespace nla3d;
using namespace nla3d::math;

dMat test_func_dmat()
{
	return dMat((uint16) 3, (uint16) 3, (double) 1.1f, 4.4f, 7.7f, 5.5f, 2.2f, 8.8f, 6.6f, 9.9f, 3.3f);
}

bool compare_double(double d1, double d2) {
  if (fabs(d1-d2) > 1.0e-5) {
    return false;
  }
  return true;
} 

int main()
{
	Mat<3,3> mat(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f);
  //test1
	Mat<3,3> mat2 = mat;
	cout << "mat=mat2:" << endl;
	cout << mat << endl;
  if (!mat.compare(mat2)) {
    LOG(ERROR) << "mat!=mat2";
    exit(1);
  }
  //test2
	cout << "mat2.zero()" << endl;
	mat2.zero();
	cout << mat2 << endl;
  if (!mat2.compare(Mat<3,3>(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f))) {
    LOG(ERROR) << "mat2!=zeros";
    exit(1);
  }

	cout << "mat2.Identity()" << endl;
	mat2.Identity();
	//mat=mat2;
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
  //test Mat<3,3>
	Mat<3,3> mat33(1.0, 2.0, 3.0, 43.0, 5.0, 6.0, 7.0, 8.0, 20.0);
	cout << "mat33:" << endl << mat33 << endl;
	cout << "det(mat33) = " << mat33.det() << "  (right result = -657)" << endl; 
  if (!compare_double(mat33.det(), -657.0f)) {
    LOG(ERROR) << "mat33.det != -657.0";
    exit(1);
  }
	cout << "inv(mat33) = " << endl << mat33.inv(mat33.det()) << endl;
  Mat<3,3> inv_mat33(-0.07914764,  0.02435312,  0.00456621,
        1.24505327,  0.00152207, -0.18721461, -0.47031963, -0.00913242,  0.12328767);
  cout << "answer inv_mat33 = " << endl << inv_mat33 << endl;
  if (!(mat33.inv(mat33.det())).compare(inv_mat33)) {
    LOG(ERROR) << "mat33.inv() is not correct";
    exit(1);
  } 

  // test Mat<4,4>
	Mat<4,4> mat44(1.0, 2.0, 3.0, 4.0, -5.0, 6.0, 7.0, 80.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 20.0);
	cout << "mat44:" << endl << mat44 << endl;
	cout << "det(mat44) = " << mat44.det() << "  (right result = -320)" << endl; 
  if (!compare_double(mat44.det(), -320.0f)) {
    LOG(ERROR) << "mat44.det != -320.0";
    exit(1);
  }
	cout << "inv(mat4) = " << endl << mat44.inv(mat44.det()) << endl;
  Mat<4,4> inv_mat44(9.50000000e-01, -1.00000000e-01, -2.65000000e+00, 1.80000000e+00,
       -3.15000000e+00, 2.00000000e-01, 5.30000000e+00, -3.35000000e+00,
       1.95000000e+00, -1.00000000e-01, -2.15000000e+00, 1.30000000e+00,
       1.25000000e-01,   6.89552582e-17,  -3.75000000e-01, 2.50000000e-01);
  cout << "answer inv_mat44 = " << endl << inv_mat44 << endl;
  if (!(mat44.inv(mat44.det())).compare(inv_mat44)) {
    LOG(ERROR) << "mat44.inv() is not correct";
    exit(1);
  } 

  // test Mat<5,5>
	Mat<5,5> mat55(1.0, 1.0, 3.0, 4.0, 5.0, 6.0, 7.0, 81.0, 9.0, 12.0, 11.0, 12.0,
      13.0, 14.0, 15.0, 16.0, 3.0, 18.0, 19.0, 20.0, 19.0, 22.0, 23.0, 24.0, 26.0);
	cout << "mat55:" << endl << mat55 << endl;
	cout << "det(mat55) = " << mat55.det() << "  (right result = 1.001E4)" << endl; 
  if (!compare_double(mat55.det(), 10005.0f)) {
    LOG(ERROR) << "mat55.det != 10005.0";
    exit(1);
  }
	cout << "inv(mat55) = " << endl << mat55.inv(mat55.det()) << endl;
  Mat<5,5> inv_mat55(-3.00649675e-01,  -1.44927536e-02,  -2.00239880e+00,
          1.72613693e-01,   1.08695652e+00,  -3.44827586e-02,
         -6.03330945e-18,   1.03448276e-01,  -6.89655172e-02,
          5.36294173e-16,  -1.65917041e-02,   1.44927536e-02,
          1.54122939e-01,  -7.09645177e-03,  -8.69565217e-02,
         -1.09045477e-01,   2.89855072e-02,   7.73583208e+00,
         -3.65917041e-01,  -4.17391304e+00,   3.64217891e-01,
         -2.89855072e-02,  -5.90134933e+00,   2.76261869e-01,
          3.17391304e+00);
  if (!(mat55.inv(mat55.det())).compare(inv_mat55)) {
    LOG(ERROR) << "mat55.inv() is not correct";
    exit(1);
  } 

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


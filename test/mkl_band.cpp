#include "sys.h"
#include "math\Mat.h"
#include "math\Mat_Band.h"
#include <mkl.h>

bool compare_double(double d1, double d2) {
  if (fabs(d1-d2) > 1.0e-6) {
    return false;
  }
  return true;
} 

void main() {
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
  Vec<5> res(9.55252e-005, -0.139571, 0.281898, 0.345268, 0.0561808);
	//второй метод
	dsptrf("U", &n, mat.Ptr(), ipv, &info);
	dsptrs("U", &n, &nrhs,  mat.Ptr(), ipv, B.Ptr(), &ldb, &info);
	cout << B << endl;
  for (uint16 i = 1; i < 6; i++) {
    if (!compare_double(B[i], res[i-1])) {
      error("Results of dsptrf is incorrect");
    }
  }
	delete[] ipv;
}

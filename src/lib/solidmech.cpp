// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "solidmech.h"

namespace nla3d {
namespace solidmech {

double J_C (const double* C) {
	return  sqrt( C[M_XX]*(C[M_YY]*C[M_ZZ]-C[M_YZ]*C[M_YZ])
              - C[M_XY]*(C[M_XY]*C[M_ZZ]-C[M_YZ]*C[M_XZ])
              + C[M_XZ]*(C[M_XY]*C[M_YZ]-C[M_YY]*C[M_XZ]) );
}

void invC_C (const double* C, const double J, double* invC) {
	double oo = 1.0/(J*J);
	invC[M_XX] = oo*(C[M_YY]*C[M_ZZ]-C[M_YZ]*C[M_YZ]);
	invC[M_XY] = oo*(C[M_XZ]*C[M_YZ]-C[M_XY]*C[M_ZZ]);
	invC[M_XZ] = oo*(C[M_XY]*C[M_YZ]-C[M_XZ]*C[M_YY]); 
	invC[M_YY] = oo*(C[M_XX]*C[M_ZZ]-C[M_XZ]*C[M_XZ]);
	invC[M_YZ] = oo*(C[M_XY]*C[M_XZ]-C[M_XX]*C[M_YZ]);
	invC[M_ZZ] = oo*(C[M_XX]*C[M_YY]-C[M_XY]*C[M_XY]);
}


void E_C (const double* C, double* E) {
  E[M_XX] += (C[M_XX]-1.0)*0.5;
  E[M_XY] += C[M_XY]*0.5;
  E[M_XZ] += C[M_XZ]*0.5;
  E[M_YY] += (C[M_YY]-1.0)*0.5;
  E[M_YZ] += C[M_YZ]*0.5;
  E[M_ZZ] += (C[M_ZZ]-1.0)*0.5;
}


void IC_C (const double* C, double* IC) {
	IC[0] = C[M_XX]+C[M_YY]+C[M_ZZ]; 
  //IC2 with star = 0.5 * IC[0]^2 - C:C
	IC[1] = C[M_XX]*C[M_YY] + C[M_YY]*C[M_ZZ] + C[M_XX]*C[M_ZZ] - C[M_XY]*C[M_XY] - C[M_YZ]*C[M_YZ] - C[M_XZ]*C[M_XZ];
	IC[2] = J_C(C);
}

} // namespace nla3d
} // namespace solidmech

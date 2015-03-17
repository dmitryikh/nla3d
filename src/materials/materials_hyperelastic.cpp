// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "materials/materials_hyperelastic.h"

namespace nla3d {
using namespace solidmech;
const double Mat_Hyper_Isotrop_General::II[6][6] = {	{1.0,0.0,0.0,0.0,0.0,0.0},
														{0.0,0.5,0.0,0.0,0.0,0.0},
														{0.0,0.0,0.5,0.0,0.0,0.0},
														{0.0,0.0,0.0,1.0,0.0,0.0},
														{0.0,0.0,0.0,0.0,0.5,0.0},
														{0.0,0.0,0.0,0.0,0.0,1.0} };

//---------------------------------------------------------
//---------------Mat_Hyper_Isotrop_General-----------------
//---------------------------------------------------------
void Mat_Hyper_Isotrop_General::getS_U (uint16 ncomp, const  solidmech::tensorComponents* comps, const double* C, double *S) {
	error("not now");
}		


void Mat_Hyper_Isotrop_General::getD_U (uint16 ncomp, const  solidmech::tensorComponents* comps, const double* C, double *D) {
	error("not now");
}	
	

//C - all times we must have 6 components
void Mat_Hyper_Isotrop_General::getS_UP (uint16 ncomp, const  solidmech::tensorComponents* comps, const double* C, double *S) {
	double press = GlobStates.getdouble(States::DOUBLE_HYDPRES);
	double alpha[2];

  double IC[3];
  solidmech::IC_C(C, IC);

	double _13I1C = 1.0/3.0*IC[0];
	double _23I2C = 2.0/3.0*IC[1];
	
	double J = IC[2];
	
	double oo = 1.0/(J*J);
	double pp = pow(J,-2.0/3.0);
	double pppp = pp*pp;
	
	double C_inv[6];
  solidmech::invC_C(C, J, C_inv);
		
	W_first_derivatives(IC[0]*pp, IC[1]*pppp, 0.0, alpha); //only first derivatives
	double ko1 = 2*alpha[AL_1]*pp;
	double ko2 = 2*alpha[AL_2]*pppp;
	
	double A[] = {I[0]-_13I1C*C_inv[0], I[1]-_13I1C*C_inv[1], I[2]-_13I1C*C_inv[2], 
					I[3]-_13I1C*C_inv[3], I[4]-_13I1C*C_inv[4], I[5]-_13I1C*C_inv[5]};
	double B[] = {IC[0]*I[0]-C[0]-_23I2C*C_inv[0], IC[0]*I[1]-C[1]-_23I2C*C_inv[1], IC[0]*I[2]-C[2]-_23I2C*C_inv[2],
					IC[0]*I[3]-C[3]-_23I2C*C_inv[3], IC[0]*I[4]-C[4]-_23I2C*C_inv[4], IC[0]*I[5]-C[5]-_23I2C*C_inv[5]};
	solidmech::tensorComponents ij;
	for (uint16 i=0; i < ncomp; i++)
	{
		ij = comps[i];
		S[i] = ko1*A[ij]+ko2*B[ij]+press*J*C_inv[ij];
	}		
}


void Mat_Hyper_Isotrop_General::getDdDp_UP (uint16 ncomp, const  solidmech::tensorComponents* comps, const double* C, double *Dd, double *Dp)
{
	double press = GlobStates.getdouble(States::DOUBLE_HYDPRES);
	
	double alpha[5];
  double IC[3];
  solidmech::IC_C(C, IC);

	double _13I1C = 1.0/3.0*IC[0];
	double _23I2C = 2.0/3.0*IC[1];
	
	double J = IC[2];
	
	double pp = pow(J,-2.0/3.0);
	double oo = pow(J,-2);
	double pppp = pp*pp;
	
	double C_inv[6];
  solidmech::invC_C(C, J, C_inv);

	/*
	AL_1	=	0
	AL_2	=	1
	AL_11	=	2
	AL_12	=	3
	AL_22	=	4
	*/
	W_second_derivatives(IC[0]*pp, IC[1]*pppp, 1.0, alpha);
	
	double A[] = {I[0]-_13I1C*C_inv[0], I[1]-_13I1C*C_inv[1], I[2]-_13I1C*C_inv[2], 
					I[3]-_13I1C*C_inv[3], I[4]-_13I1C*C_inv[4], I[5]-_13I1C*C_inv[5]};
	double B[] = {IC[0]*I[0]-C[0]-_23I2C*C_inv[0], IC[0]*I[1]-C[1]-_23I2C*C_inv[1], IC[0]*I[2]-C[2]-_23I2C*C_inv[2],
					IC[0]*I[3]-C[3]-_23I2C*C_inv[3], IC[0]*I[4]-C[4]-_23I2C*C_inv[4], IC[0]*I[5]-C[5]-_23I2C*C_inv[5]};
					
	double IIt[6][6] = {{-C_inv[M_XX]*C_inv[M_XX], -C_inv[M_XX]*C_inv[M_XY], -C_inv[M_XX]*C_inv[M_XZ], -C_inv[M_XY]*C_inv[M_XY], -C_inv[M_XY]*C_inv[M_XZ], -C_inv[M_XZ]*C_inv[M_XZ]},
		// 1211 = 11*12, 1212 = 0.5*(11*22 + 12*12), 1213 = 0.5*(11*23+13*12), 1222 = 12*22, 1223 = 0.5*(12*23+13*22), 1233 = 13*23
		 {-C_inv[M_XX]*C_inv[M_XY], -0.5*(C_inv[M_XX]*C_inv[M_YY]+C_inv[M_XY]*C_inv[M_XY]), -0.5*(C_inv[M_XX]*C_inv[M_YZ]+C_inv[M_XZ]*C_inv[M_XY]), -C_inv[M_XY]*C_inv[M_YY], -0.5*(C_inv[M_XY]*C_inv[M_YZ]+C_inv[M_XZ]*C_inv[M_YY]), -C_inv[M_XZ]*C_inv[M_YZ]},
		// 1311 = 11*13, 1312 = 0.5*(11*23+12*13), 1313 = 0.5*(11*33+13*13), 1322 = 12*23, 1323 = 0.5*(12*33+13*23), 1333 = 13*33
		 {-C_inv[M_XX]*C_inv[M_XZ], -0.5*(C_inv[M_XX]*C_inv[M_YZ]+C_inv[M_XY]*C_inv[M_XZ]), -0.5*(C_inv[M_XX]*C_inv[M_ZZ]+C_inv[M_XZ]*C_inv[M_XZ]), -C_inv[M_XY]*C_inv[M_YZ],  -0.5*(C_inv[M_XY]*C_inv[M_ZZ]+C_inv[M_XZ]*C_inv[M_YZ]), -C_inv[M_XZ]*C_inv[M_ZZ]},
		// 2211 = 12*12, 2212 = 12*22, 2213 = 12*23, 2222 = 22*22, 2223 = 22*23, 2233 = 23*23
		 {-C_inv[M_XY]*C_inv[M_XY], -C_inv[M_XY]*C_inv[M_YY], -C_inv[M_XY]*C_inv[M_YZ], -C_inv[M_YY]*C_inv[M_YY], -C_inv[M_YY]*C_inv[M_YZ], -C_inv[M_YZ]*C_inv[M_YZ]},
		// 2311 = 12*13, 2312 = 0.5*(12*23+22*13), 2313 = 0.5*(12*33+23*13), 2322 = 22*23, 2323 = 0.5*(22*33+23*23), 2333 = 23*33
		 {-C_inv[M_XY]*C_inv[M_XZ], -0.5*(C_inv[M_XY]*C_inv[M_YZ]+C_inv[M_YY]*C_inv[M_XZ]), -0.5*(C_inv[M_XY]*C_inv[M_ZZ]+C_inv[M_YZ]*C_inv[M_XZ]), -C_inv[M_YY]*C_inv[M_YZ], -0.5*(C_inv[M_YY]*C_inv[M_ZZ]+C_inv[M_YZ]*C_inv[M_YZ]), -C_inv[M_YZ]*C_inv[M_ZZ]},
		// 3311 = 13*13, 3312 = 13*23, 3313 = 13*33, 3322 = 23*23, 3323 = 23*33, 3333 = 33*33
		 {-C_inv[M_XZ]*C_inv[M_XZ], -C_inv[M_XZ]*C_inv[M_YZ], -C_inv[M_XZ]*C_inv[M_ZZ], -C_inv[M_YZ]*C_inv[M_YZ], -C_inv[M_YZ]*C_inv[M_ZZ], -C_inv[M_ZZ]*C_inv[M_ZZ]}
		};
					
	uint16 ind = 0;
	solidmech::tensorComponents ij;
	solidmech::tensorComponents kl;
	for (uint16 i=0; i < ncomp; i++)
	{
		ij = comps[i];
		for (uint16 j=i; j < ncomp; j++)
		{
			
			kl = comps[j];
			Dd[ind]=4.0*alpha[AL_11]*pppp*A[ij]*A[kl]+4.0*alpha[AL_12]*oo*(B[ij]*A[kl]+A[ij]*B[kl])+4.0*alpha[AL_22]*pppp*pppp*B[ij]*B[kl] - 
						4.0/3.0*alpha[AL_1]*pp*(C_inv[ij]*A[kl]+A[ij]*C_inv[kl]+_13I1C*C_inv[ij]*C_inv[kl]+IC[0]*IIt[ij][kl]) -
						8.0/3.0*alpha[AL_2]*pppp*(C_inv[ij]*B[kl]+B[ij]*C_inv[kl]+_23I2C*C_inv[ij]*C_inv[kl]-3.0/2.0*I[ij]*I[kl]+3.0/2.0*II[ij][kl]+IC[1]*IIt[ij][kl])+
						press*J*(C_inv[ij]*C_inv[kl]+2*IIt[ij][kl]);
			Dd[ind]=Dd[ind]/2.0;
			ind++;
		}
		Dp[i] = J*C_inv[ij];
	}
}


//---------------------------------------------------------
//------------------Mat_Comp_Neo_Hookean-------------------
//---------------------------------------------------------
void Mat_Comp_Neo_Hookean::W_first_derivatives (double I1, double I2, double I3, double* alpha)
{
	alpha[AL_1] = 0.5*MC[C_G];
	alpha[AL_2] = 0.0;	
}

void Mat_Comp_Neo_Hookean::W_second_derivatives (double I1, double I2, double I3, double* alpha)
{
	alpha[AL_1] = 0.5*MC[C_G];
	alpha[AL_2] = 0.0;	
	alpha[AL_11] = 0.0;	
	alpha[AL_22] = 0.0;	
	alpha[AL_12] = 0.0;	
}

double Mat_Comp_Neo_Hookean::W (double I1, double I2, double I3) {
    return 0.5*MC[C_G]*(I1-3.0);
}

double Mat_Comp_Neo_Hookean::getK() {
  return MC[C_K];
}

//---------------------------------------------------------
//-------------------Mat_Comp_Biderman---------------------
//---------------------------------------------------------
void Mat_Comp_Biderman::W_first_derivatives (double I1, double I2, double I3, double* alpha) {
	alpha[AL_1] = MC[C_C10]+2*MC[C_C20]*(I1-3.0)+3*MC[C_C30]*(I1-3.0)*(I1-3.0);
	alpha[AL_2] = MC[C_C01];	
}

void Mat_Comp_Biderman::W_second_derivatives (double I1, double I2, double I3, double* alpha) {
	alpha[AL_1] = MC[C_C10]+2*MC[C_C20]*(I1-3.0)+3*MC[C_C30]*(I1-3.0)*(I1-3.0);
	alpha[AL_2] = MC[C_C01];
	alpha[AL_11] = 2*MC[C_C20]+6*MC[C_C30]*(I1-3.0);	
	alpha[AL_22] = 0.0;	
	alpha[AL_12] = 0.0;	
}
double Mat_Comp_Biderman::W (double I1, double I2, double I3) {
    return MC[C_C10]*(I1-3.0) + MC[C_C20]*(I1-3.0)*(I1-3.0) + MC[C_C30]*(I1-3.0)*(I1-3.0)*(I1-3.0) + MC[C_C01]*(I2-3.0);
}

double Mat_Comp_Biderman::getK() {
  return MC[C_K];
}

//---------------------------------------------------------
//-----------------Mat_Comp_MooneyRivlin-------------------
//---------------------------------------------------------
void Mat_Comp_MooneyRivlin::W_first_derivatives (double I1, double I2, double I3, double* alpha)
{
	alpha[AL_1] = MC[C_C10];
	alpha[AL_2] = MC[C_C01];	
}

void Mat_Comp_MooneyRivlin::W_second_derivatives (double I1, double I2, double I3, double* alpha)
{
	alpha[AL_1] = MC[C_C10];
	alpha[AL_2] = MC[C_C01];
	alpha[AL_11] = 0.0;
	alpha[AL_22] = 0.0;	
	alpha[AL_12] = 0.0;	
}
double Mat_Comp_MooneyRivlin::W (double I1, double I2, double I3) {
    return MC[C_C10]*(I1-3.0) + MC[C_C01]*(I2-3.0);
}

double Mat_Comp_MooneyRivlin::getK() {
  return MC[C_K];
}

} // namespace nla3d

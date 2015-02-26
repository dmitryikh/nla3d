#include "Material.h"

const double Material::I[6] = {1.0,0.0,0.0,1.0,0.0,1.0};
const double Mat_Hyper_Isotrop_General::II[6][6] = {	{1.0,0.0,0.0,0.0,0.0,0.0},
														{0.0,0.5,0.0,0.0,0.0,0.0},
														{0.0,0.0,0.5,0.0,0.0,0.0},
														{0.0,0.0,0.0,1.0,0.0,0.0},
														{0.0,0.0,0.0,0.0,0.5,0.0},
														{0.0,0.0,0.0,0.0,0.0,1.0} };

const mat_comp MatCompsGlobal[6] = {M_XX, M_XY, M_XZ, M_YY, M_YZ, M_ZZ};
//---------------------------------------------------------
//----------------MATERIAL ABSTRACT CLASS------------------
//---------------------------------------------------------

Material::Material (uint16 num_c, mat_model _code)
{
	numC = num_c; 
	code = _code;
	MC = new double[num_c];
}

string Material::toString()
{
	//TODO: implement
	//string str;
	//str << getName() << ": ";
	//for (uint16 i=0; i < getNumC(); i++)
	//{
	//	str << "C[" << i+1 << "]=" << C[i]
	//	if (i+1 < numC)
	//		str << ", ";
	//}
	return string();
}


void Material::read_from_stream (istream &str)
{
	for (uint16 i=0; i < getNumC(); i++)
		str >> MC[i]; //TODO: check on errors
}

double& Material::getCstr (char* mname) {
	error("nor now");
	double dummy;
	return dummy;
}

//always 6 components!
double Material::getJ(const double* C) {
	double J = sqrt(C[M_XX]*(C[M_YY]*C[M_ZZ]-C[M_YZ]*C[M_YZ])-C[M_XY]*(C[M_XY]*C[M_ZZ]-C[M_YZ]*C[M_XZ])+C[M_XZ]*(C[M_XY]*C[M_YZ]-C[M_YY]*C[M_XZ]));
	return J;
}

void Material::getC_inv(const double* C, const double J, double* C_inv) 
{
	double oo = 1.0/(J*J);
	C_inv[M_XX] = oo*(C[M_YY]*C[M_ZZ]-C[M_YZ]*C[M_YZ]);
	C_inv[M_XY] = oo*(C[M_XZ]*C[M_YZ]-C[M_XY]*C[M_ZZ]);
	C_inv[M_XZ] = oo*(C[M_XY]*C[M_YZ]-C[M_XZ]*C[M_YY]); 
	C_inv[M_YY] = oo*(C[M_XX]*C[M_ZZ]-C[M_XZ]*C[M_XZ]);
	C_inv[M_YZ] = oo*(C[M_XY]*C[M_XZ]-C[M_XX]*C[M_YZ]);
	C_inv[M_ZZ] = oo*(C[M_XX]*C[M_YY]-C[M_XY]*C[M_XY]);
}

void Material::register_mat_const(uint16 num, ...) {
	numC = num;
	va_list vlist;
	va_start(vlist, num);
	for (uint16 i=0; i < num; i++) {
		MC_names.push_back(va_arg(vlist,char*));
	}
	MC = new double[num];
}



//---------------------------------------------------------
//---------------Mat_Hyper_Isotrop_General-----------------
//---------------------------------------------------------
void Mat_Hyper_Isotrop_General::getS_U (uint16 ncomp, const  mat_comp* comps, const double* C, double *S)
{
	error("not now");
}		
void Mat_Hyper_Isotrop_General::getD_U (uint16 ncomp, const  mat_comp* comps, const double* C, double *D) //for non linear U elements
{
	error("not now");
}	
	
//C - all times we must have 6 components
void Mat_Hyper_Isotrop_General::getS_UP (uint16 ncomp, const  mat_comp* comps, const double* C, double *S)
{
	assert(W_derivatives1);
	assert(W_derivatives2);
	
	// M_XX   0
	// M_XY   1
	// M_XZ   2
	// M_YY   3
	// M_YZ   4
	// M_ZZ   5
	double press = GlobStates.getdouble(States::DOUBLE_HYDPRES);
	double alpha[2];
	double I1C = C[M_XX]+C[M_YY]+C[M_ZZ]; 
	double I2Cz= C[M_XX]*C[M_YY] + C[M_YY]*C[M_ZZ] + C[M_XX]*C[M_ZZ] - C[M_XY]*C[M_XY] - C[M_YZ]*C[M_YZ] - C[M_XZ]*C[M_XZ];

	double _13I1C = 1.0/3.0*I1C;
	double _23I2C = 2.0/3.0*I2Cz;
	double J = getJ(C);
	
	double oo = 1.0/(J*J);
	double pp = pow(J,-2.0/3.0);
	double pppp = pp*pp;
	
	double C_inv[6];
	getC_inv(C,J,C_inv);
		
	W_derivatives1(I1C*pp, I2Cz*pppp, 0.0, MC, alpha); //only first derivatives
	double ko1 = 2*alpha[AL_1]*pp;
	double ko2 = 2*alpha[AL_2]*pppp;
	
	double A[] = {I[0]-_13I1C*C_inv[0], I[1]-_13I1C*C_inv[1], I[2]-_13I1C*C_inv[2], 
					I[3]-_13I1C*C_inv[3], I[4]-_13I1C*C_inv[4], I[5]-_13I1C*C_inv[5]};
	double B[] = {I1C*I[0]-C[0]-_23I2C*C_inv[0], I1C*I[1]-C[1]-_23I2C*C_inv[1], I1C*I[2]-C[2]-_23I2C*C_inv[2],
					I1C*I[3]-C[3]-_23I2C*C_inv[3], I1C*I[4]-C[4]-_23I2C*C_inv[4], I1C*I[5]-C[5]-_23I2C*C_inv[5]};
	mat_comp ij;
	for (uint16 i=0; i < ncomp; i++)
	{
		ij = comps[i];
		S[i] = ko1*A[ij]+ko2*B[ij]+press*J*C_inv[ij];
	}		
}

void Mat_Hyper_Isotrop_General::getDdDp_UP (uint16 ncomp, const  mat_comp* comps, const double* C, double *Dd, double *Dp)
{
	assert(W_derivatives1);
	assert(W_derivatives2);
	double press = GlobStates.getdouble(States::DOUBLE_HYDPRES);
	
	double alpha[5];
	double I1C = C[M_XX]+C[M_YY]+C[M_ZZ]; 
	double I2Cz= C[M_XX]*C[M_YY] + C[M_YY]*C[M_ZZ] + C[M_XX]*C[M_ZZ] - C[M_XY]*C[M_XY] - C[M_YZ]*C[M_YZ] - C[M_XZ]*C[M_XZ];
	
	double _13I1C = 1.0/3.0*I1C;
	double _23I2C = 2.0/3.0*I2Cz;
	
	double J = getJ(C);
	
	double pp = pow(J,-2.0/3.0);
	double oo = pow(J,-2);
	double pppp = pp*pp;
	
	double C_inv[6];/* = {oo*(C[M_YY]*C[M_ZZ]-C[M_YZ]*C[M_YZ]), oo*(C[M_XZ]*C[M_YZ]-C[M_XY]*C[M_ZZ]), oo*(C[M_XY]*C[M_YZ]-C[M_XZ]*C[M_YY]), 
						oo*(C[M_XX]*C[M_ZZ]-C[M_XZ]*C[M_XZ]), oo*(C[M_XY]*C[M_XZ]-C[M_XX]*C[M_YZ]), oo*(C[M_XX]*C[M_YY]-C[M_XY]*C[M_XY])};*/

	getC_inv(C, J, C_inv);

	/*
	AL_1	=	0
	AL_2	=	1
	AL_11	=	2
	AL_12	=	3
	AL_22	=	4
	*/
	W_derivatives2(I1C*pp, I2Cz*pppp, 0.0, MC, alpha); //first and second derivatives
	
	double A[] = {I[0]-_13I1C*C_inv[0], I[1]-_13I1C*C_inv[1], I[2]-_13I1C*C_inv[2], 
					I[3]-_13I1C*C_inv[3], I[4]-_13I1C*C_inv[4], I[5]-_13I1C*C_inv[5]};
	double B[] = {I1C*I[0]-C[0]-_23I2C*C_inv[0], I1C*I[1]-C[1]-_23I2C*C_inv[1], I1C*I[2]-C[2]-_23I2C*C_inv[2],
					I1C*I[3]-C[3]-_23I2C*C_inv[3], I1C*I[4]-C[4]-_23I2C*C_inv[4], I1C*I[5]-C[5]-_23I2C*C_inv[5]};
					
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
	mat_comp ij;
	mat_comp kl;
	for (uint16 i=0; i < ncomp; i++)
	{
		ij = comps[i];
		for (uint16 j=i; j < ncomp; j++)
		{
			
			kl = comps[j];
			Dd[ind]=4.0*alpha[AL_11]*pppp*A[ij]*A[kl]+4.0*alpha[AL_12]*oo*(B[ij]*A[kl]+A[ij]*B[kl])+4.0*alpha[AL_22]*pppp*pppp*B[ij]*B[kl] - 
						4.0/3.0*alpha[AL_1]*pp*(C_inv[ij]*A[kl]+A[ij]*C_inv[kl]+_13I1C*C_inv[ij]*C_inv[kl]+I1C*IIt[ij][kl]) -
						8.0/3.0*alpha[AL_2]*pppp*(C_inv[ij]*B[kl]+B[ij]*C_inv[kl]+_23I2C*C_inv[ij]*C_inv[kl]-3.0/2.0*I[ij]*I[kl]+3.0/2.0*II[ij][kl]+I2Cz*IIt[ij][kl])+
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
void Mat_Comp_Neo_Hookean::W_first_derivatives (double I1, double I2, double I3, double* mat_consts, double* alpha)
{
	alpha[AL_1] = 0.5*mat_consts[C_C10];
	alpha[AL_2] = 0.0;	
}

void Mat_Comp_Neo_Hookean::W_second_derivatives (double I1, double I2, double I3, double* mat_consts, double* alpha)
{
	alpha[AL_1] = 0.5*mat_consts[C_C10];
	alpha[AL_2] = 0.0;	
	alpha[AL_11] = 0.0;	
	alpha[AL_22] = 0.0;	
	alpha[AL_12] = 0.0;	
}

double Mat_Comp_Neo_Hookean::getGxy0()
{
	//return MC[C_C10]/(2*(1+MC[C_K])); 
	return 0.5*MC[C_C10];
}

double Mat_Comp_Neo_Hookean::getK0()
{
	//return MC[C_C10]/(3*(1-2*MC[C_K]));
	return MC[C_K];
}

double Mat_Comp_Neo_Hookean::getEx0()
{
	return MC[C_C10];
}

double Mat_Comp_Neo_Hookean::getMuxy0()
{
	return MC[C_K];
}


//---------------------------------------------------------
//-------------------Mat_Comp_Biderman---------------------
//---------------------------------------------------------
void Mat_Comp_Biderman::W_first_derivatives (double I1, double I2, double I3, double* mat_consts, double* alpha)
{
	alpha[AL_1] = mat_consts[C_C10]+2*mat_consts[C_C20]*(I1-3.0)+3*mat_consts[C_C30]*(I1-3.0)*(I1-3.0);
	alpha[AL_2] = mat_consts[C_C01];	
}

void Mat_Comp_Biderman::W_second_derivatives (double I1, double I2, double I3, double* mat_consts, double* alpha)
{
	alpha[AL_1] = mat_consts[C_C10]+2*mat_consts[C_C20]*(I1-3.0)+3*mat_consts[C_C30]*(I1-3.0)*(I1-3.0);
	alpha[AL_2] = mat_consts[C_C01];
	alpha[AL_11] = 2*mat_consts[C_C20]+6*mat_consts[C_C30]*(I1-3.0);	
	alpha[AL_22] = 0.0;	
	alpha[AL_12] = 0.0;	
}

double Mat_Comp_Biderman::getGxy0()
{
	//return MC[C_C10]/(2*(1+MC[C_K])); 
	error("why im here?");
	return 0.5*MC[C_C10];
}

double Mat_Comp_Biderman::getK0()
{
	//return MC[C_C10]/(3*(1-2*MC[C_K]));
	return MC[C_K];
}

double Mat_Comp_Biderman::getEx0()
{
	error("why im here?");
	return MC[C_C10];
}

double Mat_Comp_Biderman::getMuxy0()
{
	error("why im here?");
	return MC[C_K];
}


//---------------------------------------------------------
//-----------------Mat_Comp_MooneyRivlin-------------------
//---------------------------------------------------------
void Mat_Comp_MooneyRivlin::W_first_derivatives (double I1, double I2, double I3, double* mat_consts, double* alpha)
{
	alpha[AL_1] = mat_consts[C_C10];
	alpha[AL_2] = mat_consts[C_C01];	
}

void Mat_Comp_MooneyRivlin::W_second_derivatives (double I1, double I2, double I3, double* mat_consts, double* alpha)
{
	alpha[AL_1] = mat_consts[C_C10];
	alpha[AL_2] = mat_consts[C_C01];
	alpha[AL_11] = 0.0;
	alpha[AL_22] = 0.0;	
	alpha[AL_12] = 0.0;	
}

double Mat_Comp_MooneyRivlin::getGxy0()
{
	//return MC[C_C10]/(2*(1+MC[C_K])); 
	error("why im here?");
	return 0.5*MC[C_C10]; //TODO: be careful! this answer is wrong!
}

double Mat_Comp_MooneyRivlin::getK0()
{
	return MC[C_K];
}

double Mat_Comp_MooneyRivlin::getEx0()
{
	error("why im here?");
	return MC[C_C10];
}

double Mat_Comp_MooneyRivlin::getMuxy0()
{
	error("why im here?");
	return MC[C_K];
}


Material::matId Material::matName2matId (string matName) {
  for (uint16 i = 0; i < Material::LAST; i++) {
    if (matName.compare(Material::matModelLabels[i]) == 0) {
      return i;
    }
  }
  return Material::NOT_DEFINED;
}

Material* Material::createMaterial (string matName) {
  uint16 matId = matName2matId(matName);
  Material* mat;
  if (matId == Material::NOT_DEFINED)
    error("createMaterial: can't find material %s", matName.c_str());
  switch (matId) {
    case Material::NEO_HOOKEAN_COMP:
      mat = new Mat_Comp_Neo_Hookean();
      break;
    case Material::BIDERMAN_COMP:
      mat = new Mat_Comp_Biderman();
      break;
    case Material::MOONEYRIVLIN_COMP:
      mat = new Mat_Comp_MooneyRivlin();
      break;
    default:
      error("createMaterial: don't have a material with id %d", matId);
  }
  return mat;
}

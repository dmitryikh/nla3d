#include "Material.h"


//---------------------------------------------------------
//----------------MATERIAL ABSTRACT CLASS------------------
//---------------------------------------------------------

Material::Material (uint16 num_c, mat_model _code)
{
	numC = num_c; 
	code = _code;
	C = new double[num_c];
}

dMat Material::mat_D_d(uint16 an_type) 
{
	error("Material %s doesn't support mat_D_d creation", mat_model_labels[code]);
	return dMat(0,0);
}

dMat Material::mat_E_c(uint16 an_type, const double* C, double p_e)
{
	error("Material %s doesn't support mat_E_c creation", mat_model_labels[code]);
	return dMat(0,0);
}

dMat Material::mat_E_p(uint16 an_type, const double* C) 
{
	error("Material %s doesn't support mat_E_p creation", name.c_str());
	return dMat(0,0);
}

dMat Material::getS(uint16 an_type, const double* C, double p_e)
{
	error("Material %s doesn't support getS", name.c_str());
	return dMat(0,0);
}

string Material::toString()
{
	/*string str;
	str << getName() << ": ";
	for (uint16 i=0; i < getNumC(); i++)
	{
		str << "C[" << i+1 << "]=" << C[i]
	  if (i+1 < numC)
			str << ", ";
	}*/
	return string();
}


void Material::read_from_stream (istream &str)
{
	for (uint16 i=0; i < getNumC(); i++)
		str >> C[i]; //TODO: check on errors
}


double Material::getJ(uint16 an_type, const double* C)
{
	double J;
	switch(an_type)
	{
	case ANALYSIS_3D:
		J = sqrt(Mat<3,3>(C[0],C[3],C[5],
				  C[3],C[1],C[4],
				  C[5],C[4],C[2]).det());
		break;
	case ANALYSIS_2D_PLANE_STRAIN:
		J = sqrt(C[0]*C[1]-C[2]*C[2]);
		break;
	default:
		warning("Material::getJ analysis type %d not supported", an_type);
		J = 0.0;
	}
	return J;
}

//---------------------------------------------------------
//----------------HOOKEAN MATERIAL CLASS------------------
//---------------------------------------------------------

double Material_Hookean::getGxy0()
{
	return C[0]/(2*(1+C[1]));  //G=E/(2*(1+mu))
}

double Material_Hookean::getK0()
{
	return C[0]/(3*(1-2*C[1])); //k=E/(3*(1-2*mu))
}

double Material_Hookean::getEx0()
{
	return C[0]; //E
}

double Material_Hookean::getMuxy0()
{
	return C[1]; //mu
}

dMat Material_Hookean::mat_D_d(uint16 an_type)
{
	double G = getGxy0();
	dMat D_d(0,0);
	if (an_type == ANALYSIS_3D)
	{
		D_d.resize(6,6);
		D_d.fill(2.0*G,0.0,0.0,0.0,0.0,0.0,
				 0.0,2.0*G,0.0,0.0,0.0,0.0,
				 0.0,0.0,2.0*G,0.0,0.0,0.0,
				 0.0,0.0,0.0,G,0.0,0.0,
				 0.0,0.0,0.0,0.0,G,0.0,
				 0.0,0.0,0.0,0.0,0.0,G);
	}
	else
		warning("Material_Nookean::mat_D_d analysis type %d not supported", an_type);
	return D_d;
}


//---------------------------------------------------------
//--------COMPRESSIBLE NEO-HOOKEAN MATERIAL CLASS----------
//---------------------------------------------------------

double Material_Comp_Neo_Hookean::getGxy0()
{
	return C[0]/(2*(1+C[1])); 
}

double Material_Comp_Neo_Hookean::getK0()
{
	return C[0]/(3*(1-2*C[1]));
}

double Material_Comp_Neo_Hookean::getEx0()
{
	return C[0];
}

double Material_Comp_Neo_Hookean::getMuxy0()
{
	return C[1];
}

dMat Material_Comp_Neo_Hookean::mat_E_c(uint16 an_type, const double* C, double p_e)
{
	double G = getGxy0();
	double J = getJ(an_type, C);
	double Ic,_43JI;
	double _13GJ = 1.0/3.0*G*pow(J, -8.0/3.0);
	double _12J = 0.5*pow(J,-3);
	Vec<6> Ci;
	dMat E_c(0,0);
	switch (an_type)
	{
	case ANALYSIS_2D_PLANE_STRAIN:
		E_c.resize(3,3);
		Ic = C[0]+C[1]+1.0f;
		_43JI = 4.0/3.0*pow(J,-2)*Ic;
		E_c[0][0] =  -_13GJ*(-_43JI*C[1]*C[1]+2*C[1])+p_e*(-_12J*C[1]*C[1]);
		E_c[0][1] =  -_13GJ*(-_43JI*C[1]*C[0]+C[0]+C[1]+Ic)+p_e*(-_12J*C[1]*C[0]+pow(J,-1));
		E_c[0][2] =  -_13GJ*(_43JI*C[1]*C[2]-C[2])+p_e*(_12J*C[1]*C[2]);
		E_c[1][1] =  -_13GJ*(-_43JI*C[0]*C[0]+2*C[0])+p_e*(-_12J*C[0]*C[0]);
		E_c[1][2] =  -_13GJ*(_43JI*C[0]*C[2]-C[2])+p_e*(_12J*C[0]*C[2]);
		E_c[2][2] =  -_13GJ*(-2*_43JI*C[2]*C[2]-Ic)+p_e*(2*_12J*C[2]*C[2]-pow(J,-1));
		E_c[1][0] =  E_c[0][1];
		E_c[2][0] =  E_c[0][2];
		E_c[2][1] =  E_c[1][2];
		break;
	case ANALYSIS_3D:
		E_c.resize(6,6);
		Ic = C[0] + C[1] + C[2];
		_43JI = 4.0/3.0*pow(J, -2)*Ic;
		Ci = Vec<6>(C[1]*C[2]-C[4]*C[4], C[0]*C[2]-C[5]*C[5], C[0]*C[1]-C[3]*C[3],
			C[5]*C[4]-C[3]*C[2], C[3]*C[5]-C[0]*C[4], C[3]*C[4]-C[1]*C[5]);
		E_c[0][0] = _13GJ*(_43JI*Ci[0]*Ci[0]-2*Ci[0]) + p_e*(-_12J*Ci[0]*Ci[0]); //E1111
		E_c[0][1] = _13GJ*(_43JI*Ci[0]*Ci[1]-Ci[0]-Ci[1]-Ic*C[2]) + p_e*(-_12J*Ci[0]*Ci[1]+C[2]/J); //E1122
		E_c[0][2] = _13GJ*(_43JI*Ci[0]*Ci[2]-Ci[0]-Ci[2]-Ic*C[1]) + p_e*(-_12J*Ci[0]*Ci[2]+C[1]/J);//E1133
		E_c[0][3] = _13GJ*(_43JI*Ci[0]*Ci[3]-Ci[3])+p_e*(-_12J*Ci[0]*Ci[3]); // 0.5*(E1112+E1121)=E1112
		E_c[0][4] = _13GJ*(_43JI*Ci[0]*Ci[4]-Ci[4]+Ic*C[4])+p_e*(-_12J*Ci[0]*Ci[4]-C[4]/J); // 0.5*(E1123+E1132)=E1123
		E_c[0][5] = _13GJ*(_43JI*Ci[0]*Ci[5]-Ci[5])+p_e*(-_12J*Ci[0]*Ci[5]);// 0.5*(E1113+E1131)=E1113

		E_c[1][1] = _13GJ*(_43JI*Ci[1]*Ci[1]-2*Ci[1]) + p_e*(-_12J*Ci[1]*Ci[1]); // E2222 
		E_c[1][2] = _13GJ*(_43JI*Ci[1]*Ci[2]-Ci[1]-Ci[2]-Ic*C[0]) + p_e*(-_12J*Ci[1]*Ci[2]+C[0]/J);; // E2233
		E_c[1][3] = _13GJ*(_43JI*Ci[3]*Ci[1]-Ci[3])+p_e*(-_12J*Ci[1]*Ci[3]);// 0.5*(E2212+E2221)=E2212
		E_c[1][4] = _13GJ*(_43JI*Ci[4]*Ci[1]-Ci[4])+p_e*(-_12J*Ci[1]*Ci[4]);// 0.5*(E2223+E2232)=E2223
		E_c[1][5] = _13GJ*(_43JI*Ci[5]*Ci[1]-Ci[5]+Ic*C[5])+p_e*(-_12J*Ci[1]*Ci[5]-C[5]/J);// 0.5*(E2213+E2231)=E2213

		E_c[2][2] = _13GJ*(_43JI*Ci[2]*Ci[2]-2*Ci[2]) + p_e*(-_12J*Ci[2]*Ci[2]); //E3333
		E_c[2][3] = _13GJ*(_43JI*Ci[3]*Ci[2]-Ci[3]+Ic*C[3])+p_e*(-_12J*Ci[2]*Ci[3]-C[3]/J);// 0.5*(E3312+E3321)=E3312
		E_c[2][4] = _13GJ*(_43JI*Ci[4]*Ci[2]-Ci[4])+p_e*(-_12J*Ci[2]*Ci[4]);// 0.5*(E3323+E3332)=E3323
		E_c[2][5] = _13GJ*(_43JI*Ci[5]*Ci[2]-Ci[5])+p_e*(-_12J*Ci[2]*Ci[5]);// 0.5*(E3313+E3331)=E3313

		E_c[3][3] = 0.5*_13GJ*(2.0*_43JI*Ci[3]*Ci[3]+Ic*C[2])+0.5*p_e*(-2.0*_12J*Ci[3]*Ci[3]-C[2]/J); // 0.5*(E1212+E1221)
		E_c[3][4] = 0.5*_13GJ*(2.0*_43JI*Ci[4]*Ci[3]-Ic*C[5])+0.5*p_e*(-2.0*_12J*Ci[3]*Ci[4]+C[5]/J); // 0.5*(E1223+E1232)
		E_c[3][5] = 0.5*_13GJ*(2.0*_43JI*Ci[5]*Ci[3]-Ic*C[4])+0.5*p_e*(-2.0*_12J*Ci[3]*Ci[5]+C[4]/J); // 0.5*(E1213+E1231)

		E_c[4][4] = 0.5*_13GJ*(2.0*_43JI*Ci[4]*Ci[4]+Ic*C[0])+0.5*p_e*(-2.0*_12J*Ci[4]*Ci[4]-C[0]/J); // 0.5*(E2323+E2332)
		E_c[4][5] = 0.5*_13GJ*(2.0*_43JI*Ci[4]*Ci[5]-Ic*C[3])+0.5*p_e*(-2.0*_12J*Ci[4]*Ci[5]+C[3]/J); // 0.5*(E2331+E2313)

		E_c[5][5] = 0.5*_13GJ*(2.0*_43JI*Ci[5]*Ci[5]+Ic*C[1])+0.5*p_e*(-2.0*_12J*Ci[5]*Ci[5]-C[1]/J); // 0.5*(E1313+E1331)

		//sym.
		E_c[1][0] = E_c[0][1];
		E_c[2][0] = E_c[0][2];
		E_c[2][1] = E_c[1][2];
		E_c[3][0] = E_c[0][3];
		E_c[3][1] = E_c[1][3];
		E_c[3][2] = E_c[2][3];
		E_c[4][0] = E_c[0][4];
		E_c[4][1] = E_c[1][4];
		E_c[4][2] = E_c[2][4];
		E_c[4][3] = E_c[3][4];
		E_c[5][0] = E_c[0][5];
		E_c[5][1] = E_c[1][5];
		E_c[5][2] = E_c[2][5];
		E_c[5][3] = E_c[3][5];
		E_c[5][4] = E_c[4][5];
		break;
	default:
		warning("Material_Comp_Neo_Hookean::mat_E_c analysis type %d not supported", an_type);
	}
	return E_c;
}

dMat Material_Comp_Neo_Hookean::getS(uint16 an_type, const double* C, double p_e)
{
	assert(C);
	dMat S(0,0);
	double J = getJ(an_type, C);
	double G = getGxy0();
	double GJ;
	Vec<6> Ci;
	double Ic;
	switch (an_type)
	{
	case ANALYSIS_2D_PLANE_STRAIN:
		Ic = C[0]+C[1]+1.0f;
		S.resize(3,1);
		GJ = G*pow(J,-2.0/3.0);
		//TODO: проверить
		S[0][0] = GJ*(1.0-1.0/3.0*Ic*pow(J,-2)*C[1])+p_e/J*C[1];
		S[1][0] = GJ*(1.0-1.0/3.0*Ic*pow(J,-2)*C[0])+p_e/J*C[0];
		S[2][0] = GJ*(+1.0/3.0*Ic*pow(J,-2)*C[2])-p_e/J*C[2];
		break;
	case ANALYSIS_3D:
		Ic = C[0] + C[1] + C[2];
		//вектор компонент матрицы алг. дополнений C11,C22,C33,C12,C23,C13
		Ci = Vec<6>(C[1]*C[2]-C[4]*C[4], C[0]*C[2]-C[5]*C[5], C[0]*C[1]-C[3]*C[3],
				C[5]*C[4]-C[3]*C[2], C[3]*C[5]-C[0]*C[4], C[3]*C[4]-C[1]*C[5]);
		S.resize(6,1);
		S.cpVec<6>(Ci*(-1.0/3.0*Ic*G*pow(J,-8.0/3.0)+p_e/J) + Vec<6>(1.0,1.0,1.0,0.0,0.0,0.0)*(G*pow(J,-2.0/3.0)));
		break;
	default:
		warning("Material_Comp_Neo_Hookean::getS analysis type %d not supported", an_type);
	}
	return S;
}

dMat Material_Comp_Neo_Hookean::mat_E_p(uint16 an_type, const double* C)
{
	assert(C);
	
	dMat E_p(0,0);
	double J = getJ(an_type, C);
	double inv_J = 1.0/J;
	Vec<6> Ci;

	switch (an_type)
	{
	case ANALYSIS_2D_PLANE_STRAIN:
		E_p.resize(3,1);
		E_p[0][0] =  inv_J*C[1];
		E_p[1][0] =  inv_J*C[0];
		E_p[2][0] =  -inv_J*C[2];
		break;
	case ANALYSIS_3D:
		E_p.resize(6,1);
		Ci = Vec<6>(C[1]*C[2]-C[4]*C[4], C[0]*C[2]-C[5]*C[5], C[0]*C[1]-C[3]*C[3],
				C[5]*C[4]-C[3]*C[2], C[3]*C[5]-C[0]*C[4], C[3]*C[4]-C[1]*C[5]);
		E_p.cpVec<6>(Ci*inv_J);
		break;
	default:
		warning("Material_Comp_Neo_Hookean::mat_E_p analysis type %d not supported", an_type);
	}
	return E_p;
}


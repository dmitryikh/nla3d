#include "materials/material.h"

const double Material::I[6] = {1.0,0.0,0.0,1.0,0.0,1.0};

//---------------------------------------------------------
//----------------MATERIAL ABSTRACT CLASS------------------
//---------------------------------------------------------
Material::Material (uint16 num_c)
{
  code = 0;
	numC = num_c; 
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


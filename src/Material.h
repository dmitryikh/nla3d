#pragma once
#include <string>
#include "math\Vec.h"
#include "math\Mat.h"

enum mat_model
{
	MAT_NOT_DEFINED,
	MAT_HOOKEAN,
	MAT_COMP_NEO_HOOKEAN
};
const char* const mat_model_labels[]={"UNDEFINED", "Hookean", "Compressible Neo-Hookean"};

class Material
{
public:
	Material () : C(NULL), numC(0), code(MAT_NOT_DEFINED)
	{ 
	}
	Material (uint16 num_c, mat_model _code);
	~Material()   // TODO: discover the virtual destructor
	{
		if (C) delete[] C; 

		C = NULL;
	}
	virtual double getK0 () = 0; //initial(linearized) Bulk module
	virtual double getGxy0 () = 0; //initial Shear module
	virtual double getEx0 () = 0; //initial tension/comp. module
	virtual double getMuxy0 () = 0; //initial Poisson coef.
	virtual dMat mat_D_d(uint16 an_type); // для лин. задач, связь между девиаторами напр и деф.
	virtual dMat mat_E_c(uint16 an_type, const double* C, double p_e); // для нелин. задач смешанного метода: функция возвращ. касат. матрицу упругости (дифф. по C)
	virtual dMat mat_E_p(uint16 an_type, const double* C); // для нелин. задач смешанного метода, касательный модуль по гидростатическому давлению
	virtual dMat getS(uint16 an_type, const double* C, double p_e); //получить напряжения П-К2 из компонент мер деформаций
	virtual string toString();
	string getName();
	uint16 getCode ();
	double& Ci (uint16 i);
	uint16 getNumC ();
	void read_from_stream (istream &str);
	double getJ(uint16 an_type, const double* C);
protected:
	double* C;
	uint16 numC;
	mat_model code;
	string name;
	
};

//базовый класс для гиперупругих материалов
class Material_Hyper : public Material
{
public:
	dMat mat_E(uint16 an_type, const double* C, double* par); // для нелин. задач смешанного метода: функция возвращ. касат. матрицу упругости (дифф. по C)
	dMat mat_Ep(uint16 an_type, const double* C, double* par); // для нелин. задач смешанного метода, касательный модуль по гидростатическому давлению
	dMat vec_S(uint16 an_type, const double* C, double* par); //получить напряжения П-К2 из компонент мер деформаций
	virtual void W_first_derivatives (double I1, double I2, double J, double *ksi)=0; //ksi1, ksi2, ksiJ
	virtual void W_derivatives (double I1, double I2, double J, double *ksi)=0;// ksi1, ksi2, ksi11, ksi12, ksi22, ksiJ, ksiJJ
	virtual void get_Eijkl (uint16 an_type, const double* C, double* par, double *E);
	virtual void get_S (uint16 an_type, const double* C, double* par, double *S);
	void get_matC(uint16 an_type, const double* C, Mat<3,3> &mat);
	// E1111, E1122, E1133, 0.5E1112+1121
	static const uint16 mat_index[6][2]; //матрица индексов 11 22 33 12 23 13
	static const uint16 pds_index[3];
	static const Mat<3,3> matI;
};

//базовый класс для гиперупругих материалов для смешанной формулировки, W_vol = 0.5k(J-1)^2;
class Material_Hyper_Mixed : public Material_Hyper
{
public:
	virtual void W_first_derivatives (double I1, double I2, double J, double *ksi)=0; //ksi1, ksi2
	virtual void W_derivatives (double I1_, double I2_, double *ksi)=0;// ksi1, ksi2, ksi11, ksi12, ksi22
	virtual void get_Eijkl (uint16 an_type, const double* C, double* par, double *E);
	virtual void get_S (uint16 an_type, const double* C, double* par, double *S);
};


// constants:
// C1 = E,  C2 = mu
class Material_Hookean : public Material
{
public:
	Material_Hookean() : Material(2, MAT_HOOKEAN) 
	{
	}
	double getK0 (); //initial(linearized) Bulk module
	double getGxy0 (); //initial Shear module
	double getEx0 (); //initial tension/comp. module
	double getMuxy0 (); //initial Poisson coef.
	dMat mat_D_d(uint16 an_type);
};



//constants:
// C1 = E; C2 = mu
class Material_Comp_Neo_Hookean : public Material
{
public:
	Material_Comp_Neo_Hookean() : Material(2, MAT_COMP_NEO_HOOKEAN) 
	{
	}
	double getK0 (); //initial(linearized) Bulk module
	double getGxy0 (); //initial Shear module
	double getEx0 (); //initial tension/comp. module
	double getMuxy0 (); //initial Poisson coef.

	virtual dMat mat_E_c(uint16 an_type, const double* C, double p_e); // для нелин. задач смешанного метода: функция возвращ. касат. матрицу упругости (дифф. по C)
	virtual dMat mat_E_p(uint16 an_type, const double* C); // для нелин. задач смешанного метода, касательный модуль по гидростатическому давлению
	virtual dMat getS(uint16 an_type, const double* C, double p_e); //получить напряжения П-К2 из компонент мер деформаций
};


inline string Material::getName()
{
	return mat_model_labels[code];
}

inline uint16 Material::getCode ()
{
	return code;
}

inline double& Material::Ci (uint16 i)
{
	assert(i < numC);
	return C[i];
}

inline uint16 Material::getNumC ()
{
	return numC;
}

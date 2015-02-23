#pragma once
#include <string>
#include "math\Vec.h"
#include "math\Mat.h"
#include "States.h"

enum mat_comp
{
	M_XX =  0,
	M_XY =  1,
	M_XZ =  2,
	M_YY =  3,
	M_YZ =  4,
	M_ZZ =  5
};

enum mat_model
{
	MAT_NOT_DEFINED,
  MAT_COMP_NEO_HOOKEAN,
  MAT_COMP_BIDERMAN,
  MAT_COMP_MOONEYRIVLIN,
  MAT_LAST
};
const char* const mat_model_labels[]={"UNDEFINED",
  "Neo-Hookean",
  "Biderman",
  "Mooney-Rivlin"
};

// Material must have named material constants

class Material
{
public:
	enum mat_func_deriv
	{
		AL_1	=	0,
		AL_2	=	1,
		AL_11	=	2,
		AL_12	=	3,
		AL_22	=	4
	};
	Material () : MC(NULL), numC(0), code(MAT_NOT_DEFINED)
	{ 
	}
	Material (uint16 num_c, mat_model _code);
	~Material()   // TODO: discover the virtual destructor
	{
		if (MC) delete[] MC; 
		MC = NULL;
	}
	virtual double getK0 () = 0; //initial(linearized) Bulk module
	virtual double getGxy0 () = 0; //initial Shear module
	virtual double getEx0 () = 0; //initial tension/comp. module
	virtual double getMuxy0 () = 0; //initial Poisson coef.
	virtual void getS_U (uint16 ncomp, const mat_comp* comps, const double* C, double *S) = 0;
	virtual void getS_UP (uint16 ncomp, const mat_comp* comps, const double* C, double *S) = 0; //
	virtual void getD_U (uint16 ncomp, const mat_comp* comps, const double* C, double *D) = 0; //for non linear U elements
	virtual void getDdDp_UP (uint16 ncomp, const mat_comp* comps, const double* C, double *Dd, double *Dp) = 0; //for nonlinear U-P elements
	virtual string toString();
	string getName();
	uint16 getCode ();
	double& Ci (uint16 i);
	double& getCstr (char* mname);
	uint16 getNumC ();
	void read_from_stream (istream &str);
	static double getJ(const double* C);
	static void getC_inv(const double* C, const double J, double* C_inv);
	static const double I[6];
protected:
	void register_mat_const(uint16 num, ...);
	double* MC;
	vector<string> MC_names;
	uint16 numC;
	mat_model code;
	string name;
	
};

class Mat_Hyper_Isotrop_General : public Material
{
	public:
	Mat_Hyper_Isotrop_General ()
	{
		W_derivatives1 = NULL;
		W_derivatives2 = NULL;	
	}
	
	virtual void getS_U (uint16 ncomp, const  mat_comp* comps, const double* C, double *S);
	virtual void getS_UP (uint16 ncomp, const  mat_comp* comps, const double* C, double *S); //
	virtual void getD_U (uint16 ncomp, const  mat_comp* comps, const double* C, double *D); //for non linear U elements
	virtual void getDdDp_UP (uint16 ncomp, const  mat_comp* comps, const double* C, double *Dd, double *Dp); //for nonlinear U-P elements	
	void getS_P_test (uint16 ncomp, const  mat_comp* comps, const double* C, double *S);
	void Mat_Hyper_Isotrop_General::getDdDp_UP_test (uint16 ncomp, const  mat_comp* comps, const double* C, double *Dd, double *Dp);
	void (*W_derivatives1) (double, double, double, double*, double*);
	void (*W_derivatives2) (double, double, double, double*, double*);
	
	static const double II[6][6];
};

class Mat_Comp_Neo_Hookean : public Mat_Hyper_Isotrop_General
{
	public:
	enum Mat_Comp_Neo_Hookean_CONST {
		C_C10 = 0,
		C_K
	};
	
	Mat_Comp_Neo_Hookean () {
		register_mat_const(2,"C10","K");
		W_derivatives1 = &W_first_derivatives;
		W_derivatives2 = &W_second_derivatives;
	}
	static void W_first_derivatives (double I1, double I2, double I3, double* mat_consts, double *alpha);
	static void W_second_derivatives (double I1, double I2, double I3, double* mat_consts, double *alpha);
	double getK0 (); //initial(linearized) Bulk module
	double getGxy0 (); //initial Shear module
	double getEx0 (); //initial tension/comp. module
	double getMuxy0 (); //initial Poisson coef.
};

class Mat_Comp_Biderman : public Mat_Hyper_Isotrop_General
{
	public:
	enum Mat_Comp_Biderman_CONST {
		C_C10 = 0,
		C_C20,
		C_C30,
		C_C01,
		C_K
	};
	
	Mat_Comp_Biderman () {
		register_mat_const(5,"C10","C20","C30","C01","K");
		W_derivatives1 = &W_first_derivatives;
		W_derivatives2 = &W_second_derivatives;
	}
	static void W_first_derivatives (double I1, double I2, double I3, double* mat_consts, double *alpha);
	static void W_second_derivatives (double I1, double I2, double I3, double* mat_consts, double *alpha);
	double getK0 (); //initial(linearized) Bulk module
	double getGxy0 (); //initial Shear module
	double getEx0 (); //initial tension/comp. module
	double getMuxy0 (); //initial Poisson coef.
};

class Mat_Comp_MooneyRivlin : public Mat_Hyper_Isotrop_General
{
	public:
	enum Mat_Comp_MooneyRivlin_CONST {
		C_C10 = 0,
		C_C01,
		C_K
	};
	
	Mat_Comp_MooneyRivlin () {
		register_mat_const(3,"C10","C01","K");
		W_derivatives1 = &W_first_derivatives;
		W_derivatives2 = &W_second_derivatives;
	}
	static void W_first_derivatives (double I1, double I2, double I3, double* mat_consts, double *alpha);
	static void W_second_derivatives (double I1, double I2, double I3, double* mat_consts, double *alpha);
	double getK0 (); //initial(linearized) Bulk module
	double getGxy0 (); //initial Shear module
	double getEx0 (); //initial tension/comp. module
	double getMuxy0 (); //initial Poisson coef.
};

// ---=== FUNCTIONS ===--- //
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
	return MC[i];
}

inline uint16 Material::getNumC ()
{
	return numC;
}

uint16 matName2matId (string matName); 
Material* createMaterial (string matName); 

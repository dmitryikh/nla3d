#pragma once
#include "materials/material.h"

//----------------------------------------------//
//-----------Mat_Hyper_Isotrop_General----------//
//----------------------------------------------//
class Mat_Hyper_Isotrop_General : public Material {
	public:
	//Mat_Hyper_Isotrop_General () {
	//}
	//for non linear U elements
	void getS_U (uint16 ncomp, const  tensorComponents* comps, const double* C, double *S);
	void getD_U (uint16 ncomp, const  tensorComponents* comps, const double* C, double *D); 
  //for nonlinear U-P elements	
	void getS_UP (uint16 ncomp, const  tensorComponents* comps, const double* C, double *S);
	void getDdDp_UP (uint16 ncomp, const  tensorComponents* comps, const double* C, double *Dd, double *Dp); 

	virtual void W_first_derivatives (double I1, double I2, double I3, double *alpha) = 0;
	virtual void W_second_derivatives (double I1, double I2, double I3, double *alpha) = 0;

  virtual double getK() = 0;
	
	static const double II[6][6];
};

class Mat_Comp_Neo_Hookean : public Mat_Hyper_Isotrop_General
{
	public:
	enum Mat_Comp_Neo_Hookean_CONST {
		C_G = 0,
		C_K
	};
	
	Mat_Comp_Neo_Hookean () {
		register_mat_const(2,"G","K");
    name = "Neo-Hookean";
	}

	void W_first_derivatives (double I1, double I2, double I3, double *alpha);
	void W_second_derivatives (double I1, double I2, double I3, double *alpha);

  double getK();
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
    name = "Biderman";
	}
	void W_first_derivatives (double I1, double I2, double I3, double *alpha);
	void W_second_derivatives (double I1, double I2, double I3, double *alpha);

  double getK();
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
    name = "Money-Rivlin";
	}
	void W_first_derivatives (double I1, double I2, double I3, double *alpha);
	void W_second_derivatives (double I1, double I2, double I3, double *alpha);

  double getK();
};


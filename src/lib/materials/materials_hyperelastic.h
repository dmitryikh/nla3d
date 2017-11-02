// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#pragma once
#include "materials/material.h"
#include "solidmech.h"
#include <Eigen/Dense>

namespace nla3d {

//----------------------------------------------//
//-----------Mat_Hyper_Isotrop_General----------//
//----------------------------------------------//
class Mat_Hyper_Isotrop_General : public Material {

  public:
  //Mat_Hyper_Isotrop_General () {
  //}
  // constants for derives in double array
  enum mat_func_deriv {
    AL_1  = 0,
    AL_2  = 1,
    AL_11 = 2,
    AL_12 = 3,
    AL_22 = 4
  };
  //for non linear U elements
  void getS_U (uint16 ncomp, const  solidmech::tensorComponents* comps, const double* C, double *S);
  void getD_U (uint16 ncomp, const  solidmech::tensorComponents* comps, const double* C, double *D);
  //for nonlinear U-P elements
  void getS_UP (uint16 ncomp, const  solidmech::tensorComponents* comps, const double* C, const double press, double *S);
  void getDdDp_UP (uint16 ncomp, const  solidmech::tensorComponents* comps, const double* C, const double press, double *Dd, double *Dp);


  virtual void W_first_derivatives (double I1, double I2, double I3, double *alpha) = 0;
  virtual void W_second_derivatives (double I1, double I2, double I3, double *alpha) = 0;
  virtual double W (double I1, double I2, double I3) = 0;

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
    double W (double I1, double I2, double I3);

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
    double W (double I1, double I2, double I3);

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
    double W (double I1, double I2, double I3);

    double getK();
};

} // namespace nla3d

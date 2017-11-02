// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#pragma once
#include "sys.h"

namespace nla3d {
namespace solidmech {

// labels for tensor components stored in 1-dim array
enum tensorComponents {
  M_XX =  0,
  M_XY =  1,
  M_XZ =  2,
  M_YY =  3,
  M_YZ =  4,
  M_ZZ =  5
};

// global order of tensor components in 1-dim array
const tensorComponents defaultTensorComponents[] = {M_XX, M_XY, M_XZ, M_YY, M_YZ, M_ZZ};

const int LeviCivita[3][3][3] = {
{ {0,0,0},
  {0,0,1},
  {0,-1,0} },
{ {0,0,-1},
  {0,0,0},
  {1,0,0} },
{ {0,1,0},
  {-1,0,0},
  {0,0,0} } };

const int I[3][3] = {
    {1,0,0},
    {0,1,0},
    {0,0,1}};

const char* const labelsTensorComponent[]={"XX","XY", "XZ", "YY", "YZ", "ZZ"};
// 3x3 symmetric tensor is stored in 1-dim array like order:
// C = [CXX, CXY, CXZ, CYY, CYZ, CZZ]
//
// get volume deformation from right Cauchy–Green deformation tensor C = F^T * F
double J_C(const double* C);
// get inverse tensor from right Cauchy–Green deformation tensor C = F^T * F
void invC_C(const double* C, const double J, double* C_inv);
// get Green strain tensor (Largrangian) E = 0.5* (C - I)
void E_C(const double* C, double* E);
// get 3 Invariants from right Cauchy-Green deformation tensor C = F^T * F
void IC_C(const double* C, double* IC);

} // namespace solidmech
} //namespace nla3d

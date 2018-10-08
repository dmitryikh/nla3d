// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"

namespace nla3d {

class ElementINTER3 : public ElementWEDGE {
public:
  ElementINTER3();

  void pre();

  void buildK();
  void make_D(Eigen::MatrixXd& D);
  Eigen::MatrixXd make_B(uint16 np, uint16 npj);
  Eigen::MatrixXd make_subB(uint16 np, uint16 npj);
  Eigen::MatrixXd make_T();

  void makeJacob();

  void update();

  uint16 getNNodes();
  uint16 getDim();

  uint16 nOfIntPoints();

  double intWeight(uint16 np);
  double intL1(uint16 np);
  double intL2(uint16 np);
  double intL3(uint16 np);

  double radoIntWeight(uint16 np, uint16 npj);
  double radoIntL1(uint16 np);
  double radoIntL2(uint16 np,uint16 npj);
  double radoIntL3(uint16 np,uint16 npj);

  // strength
  double kn = 0.0, ks = 0.0;

  // surface averaged strains
  math::Vec<3> strains;
  // surface averaged stress
  math::Vec<3> stress;
  // normal to triangle
  math::Vec<3> normal;

  uint16 i_int = 4; // index of integration scheme
  double det = 0.; // determinant of Jacob matrix

  //postproc procedures
  bool getVector(math::Vec<3>* vector, vectorQuery query, uint16 gp, const double scale);
  bool getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale);
};

inline uint16 ElementINTER3::getNNodes() {
  return 6;
}

inline uint16 ElementINTER3::getDim(){
  return 3;
}


// structures for convenient keeping of quadrature constants
struct TrianglePt {
  double L1;
  double L2;
  double L3;
  double W;
};

struct RadoPt {
  double AJ;
  double H;
  double AI;
  double AS;
};

// gauss quadrature for 3D trinangle from (0) to (1) in L-coordinates
// 1st order
static const TrianglePt _triangle_o1[] = {
  {1./3., 1./3., 1./3., 1.}
};

// 2nd order
static const TrianglePt _triangle_o2[] = {
  {1./2.,     0.,     1./2.,  1./3.},
  {1./2.,     1./2.,  0.,     1./3.},
  {0.,        1./2.,  1./2.,  1./3.}
};

// 3d order
static const TrianglePt _triangle_o3[] = {
  {1./3.,     1./3.,     1./3.,   -27./48.},
  {11./15.,   2./15.,    2./15.,   25./48.},
  {2./15.,    2./15.,    11./15.,  25./48.},
  {2./15.,    11./15.,   2./15.,   25./48.}
};

// 4 order
static const TrianglePt _triangle_o4[] = {
  {1./3.,     1./3.,     1./3.,    27./60.},
  {1./2.,     1./2.,     0.,       8./60.},
  {0.,        1./2.,     1./2.,    8./60.},
  {1./2.,     0.,        1./2.,    8./60.},
  {1.,        0.,        0.,       3./60.},
  {0.,        1.,        0.,       3./60.},
  {0.,        0.,        1.,       3./60.}
};

// 5 order
static const TrianglePt _triangle_o5[] = {
  {1./3.,     1./3.,     1./3.,    0.225 },
  {0.05971587,     0.47014206,     0.47014206, 0.13239415},
  {0.47014206,     0.05971587,     0.47014206, 0.13239415},
  {0.47014206,     0.47014206,     0.05971587, 0.13239415},
  {0.79742699,     0.10128651,     0.10128651, 0.12593918},
  {0.10128651,     0.79742699,     0.10128651, 0.12593918},
  {0.10128651,     0.10128651,     0.79742699, 0.12593918}
};

static const RadoPt _rado_o3[] = {
  {0.1127016654, 0.2777777778, 0.0885879595, 0.2204622112},
  {0.5,          0.4444444444, 0.4094668644, 0.3881934688},
  {0.8872983346, 0.2777777778, 0.7876594618, 0.3288443200}
};

// array of number of quadrature points in integration scheme
static const uint16 _np_triangle[] = {
  sizeof(_triangle_o1) / sizeof(TrianglePt),
  sizeof(_triangle_o2) / sizeof(TrianglePt),
  sizeof(_triangle_o3) / sizeof(TrianglePt),
  sizeof(_triangle_o4) / sizeof(TrianglePt),
  sizeof(_triangle_o5) / sizeof(TrianglePt)
};

inline uint16 ElementINTER3::nOfIntPoints(){
  //return _np_triangle[i_int];
  return 3;
}

static const TrianglePt* _table_triangle[] = {
  _triangle_o1,
  _triangle_o2,
  _triangle_o3,
  _triangle_o4,
  _triangle_o5
};

inline double ElementINTER3::intWeight(uint16 np){
  return _table_triangle[i_int][np].W * det;
}

inline double ElementINTER3::radoIntWeight(uint16 np, uint16 npj){
  return _rado_o3[np].AS*_rado_o3[npj].H*(1-radoIntL1(np))* det;
}

inline double ElementINTER3::intL1(uint16 np){
  return _table_triangle[i_int][np].L1;
}

inline double ElementINTER3::radoIntL1(uint16 np){
  return _rado_o3[np].AI;
}

inline double ElementINTER3::intL2(uint16 np){
  return _table_triangle[i_int][np].L2;
}

inline double ElementINTER3::radoIntL2(uint16 np, uint16 npj){
  return _rado_o3[npj].AJ*(1-radoIntL1(np));
}

inline double ElementINTER3::intL3(uint16 np){
  return _table_triangle[i_int][np].L3;
}

inline double ElementINTER3::radoIntL3(uint16 np, uint16 npj){
  return 1.-radoIntL1(np)-radoIntL2(np,npj);
}

} //namespace nla3d

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"

// Formulations for isoparametric FE for different shapes. Isoparametic formulation means usage of the
// same shape functions for geometry interpolation and for field (displacement, for instance)
//   interpolation.
// These derived classes specify next things:
// 1. Specify form functions
// 2. Specify integration rule (Gauss quadrature)
// 3. Specify field derivatives in local and global coordinates

namespace nla3d {

// structures for convenient keeping of quadrature constants
struct QuadPt1D {
  double r;
  double w;
};


struct QuadPt2D {
  double r;
  double s;
  double w;
};


struct QuadPt3D {
  double r;
  double s;
  double t;
  double w;
};

// gauss quadrature for 2D rect from (-1,-1) to (1,1)
// 1st order
const static QuadPt2D _quad_quad_o1[] = {
  {+0.0000000000000000, +0.0000000000000000, +4.0000000000000000}
};


// 2st order
const static QuadPt2D _quad_quad_o2[] = {
  {-0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {-0.5773502691896258, +0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, +0.5773502691896258, +1.0000000000000000}
};


// 2st order
const static QuadPt2D _quad_quad_o3[] = {
  {-0.7745966692414834, -0.7745966692414834, +0.3086419753086420},
  {+0.0000000000000000, -0.7745966692414834, +0.4938271604938271},
  {+0.7745966692414834, -0.7745966692414834, +0.3086419753086420},

  {-0.7745966692414834, +0.0000000000000000, +0.4938271604938271},
  {+0.0000000000000000, +0.0000000000000000, +0.7901234567901234},
  {+0.7745966692414834, +0.0000000000000000, +0.4938271604938271},

  {-0.7745966692414834, +0.7745966692414834, +0.3086419753086420},
  {+0.0000000000000000, +0.7745966692414834, +0.4938271604938271},
  {+0.7745966692414834, +0.7745966692414834, +0.3086419753086420}
};


// array of number of quadrature points in integration scheme
const static uint16 _np_quad[] = {
  sizeof(_quad_quad_o1) / sizeof(QuadPt2D),
  sizeof(_quad_quad_o2) / sizeof(QuadPt2D),
  sizeof(_quad_quad_o3) / sizeof(QuadPt2D)
};


const static QuadPt2D* _table_quad[] = {
  _quad_quad_o1,
  _quad_quad_o2,
  _quad_quad_o3
};


// gauss quadrature for 3D brick from (-1,-1,-1) to (1,1,1)
// 1st order
const static QuadPt3D _quad_hexahedron_o1[] = {
  {+0.0000000000000000, +0.0000000000000000, +0.0000000000000000, +8.0000000000000000}
};


// 2st order
const static QuadPt3D _quad_hexahedron_o2[] = {
  {-0.5773502691896258, -0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, -0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {-0.5773502691896258, +0.5773502691896258, -0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, +0.5773502691896258, -0.5773502691896258, +1.0000000000000000},

  {-0.5773502691896258, -0.5773502691896258, +0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, -0.5773502691896258, +0.5773502691896258, +1.0000000000000000},
  {-0.5773502691896258, +0.5773502691896258, +0.5773502691896258, +1.0000000000000000},
  {+0.5773502691896258, +0.5773502691896258, +0.5773502691896258, +1.0000000000000000}
};


// 3st order
const static QuadPt3D _quad_hexahedron_o3[] = {
  {-0.7745966692414834, -0.7745966692414834, -0.7745966692414834, +0.1714677640603567},
  {+0.0000000000000000, -0.7745966692414834, -0.7745966692414834, +0.2743484224965707},
  {+0.7745966692414834, -0.7745966692414834, -0.7745966692414834, +0.1714677640603567},

  {-0.7745966692414834, +0.0000000000000000, -0.7745966692414834, +0.2743484224965707},
  {+0.0000000000000000, +0.0000000000000000, -0.7745966692414834, +0.4389574759945130},
  {+0.7745966692414834, +0.0000000000000000, -0.7745966692414834, +0.2743484224965707},

  {-0.7745966692414834, +0.7745966692414834, -0.7745966692414834, +0.1714677640603567},
  {+0.0000000000000000, +0.7745966692414834, -0.7745966692414834, +0.2743484224965707},
  {+0.7745966692414834, +0.7745966692414834, -0.7745966692414834, +0.1714677640603567},


  {-0.7745966692414834, -0.7745966692414834, +0.0000000000000000, +0.2743484224965707},
  {+0.0000000000000000, -0.7745966692414834, +0.0000000000000000, +0.4389574759945130},
  {+0.7745966692414834, -0.7745966692414834, +0.0000000000000000, +0.2743484224965707},

  {-0.7745966692414834, +0.0000000000000000, +0.0000000000000000, +0.4389574759945130},
  {+0.0000000000000000, +0.0000000000000000, +0.0000000000000000, +0.7023319615912207},
  {+0.7745966692414834, +0.0000000000000000, +0.0000000000000000, +0.4389574759945130},

  {-0.7745966692414834, +0.7745966692414834, +0.0000000000000000, +0.2743484224965707},
  {+0.0000000000000000, +0.7745966692414834, +0.0000000000000000, +0.4389574759945130},
  {+0.7745966692414834, +0.7745966692414834, +0.0000000000000000, +0.2743484224965707},


  {-0.7745966692414834, -0.7745966692414834, +0.7745966692414834, +0.1714677640603567},
  {+0.0000000000000000, -0.7745966692414834, +0.7745966692414834, +0.2743484224965707},
  {+0.7745966692414834, -0.7745966692414834, +0.7745966692414834, +0.1714677640603567},

  {-0.7745966692414834, +0.0000000000000000, +0.7745966692414834, +0.2743484224965707},
  {+0.0000000000000000, +0.0000000000000000, +0.7745966692414834, +0.4389574759945130},
  {+0.7745966692414834, +0.0000000000000000, +0.7745966692414834, +0.2743484224965707},

  {-0.7745966692414834, +0.7745966692414834, +0.7745966692414834, +0.1714677640603567},
  {+0.0000000000000000, +0.7745966692414834, +0.7745966692414834, +0.2743484224965707},
  {+0.7745966692414834, +0.7745966692414834, +0.7745966692414834, +0.1714677640603567}
};


const static uint16 _np_hexahedron[] = {
  sizeof(_quad_hexahedron_o1) / sizeof(QuadPt3D),
  sizeof(_quad_hexahedron_o2) / sizeof(QuadPt3D),
  sizeof(_quad_hexahedron_o3) / sizeof(QuadPt3D)
};


const static QuadPt3D* _table_hexahedron[] = {
  _quad_hexahedron_o1,
  _quad_hexahedron_o2,
  _quad_hexahedron_o3
};


class ElementIsoParamQUAD : public ElementQUAD {
  public:
    std::vector<math::Mat<4, 2> > NiXj; //derivates form function / local coordinates
    std::vector<double> det;  //Jacobian

    // function to calculate all staff for isoparametric FE
    void makeJacob(); 
    double g_weight(uint16 np);
    double volume();
    void np2rst(uint16 np, double *xi); //by number of gauss point find local coordinates
    uint16 nOfIntPoints();
    math::Vec<4> formFunc(double r, double s);
    math::Mat<4, 2> formFuncDeriv(double r, double s);

  protected:
    uint16 i_int = 0; // index of integration scheme
};


class ElementIsoParamHEXAHEDRON : public ElementHEXAHEDRON {
  public:
    std::vector<math::Mat<8, 3> > NiXj; //derivates form function / local coordinates
    std::vector<double> det;  //Jacobian

    // function to calculate all staff for isoparametric FE
    void makeJacob(); 
    double g_weight(uint16 np);
    double volume();
    void np2rst(uint16 np, double *xi); //by number of gauss point find local coordinates
    uint16 nOfIntPoints();
    math::Vec<8> formFunc(double r, double s, double t);
    math::Mat<8, 3> formFuncDeriv(double r, double s, double t);

  protected:
    uint16 i_int = 0; // index of integration scheme
};


inline double ElementIsoParamQUAD::g_weight(uint16 np) {
  assert(np < _np_quad[i_int]);
  assert(det.size() != 0);

  return _table_quad[i_int][np].w * det[np];
}


inline void ElementIsoParamQUAD::np2rst (uint16 np, double *v) {
  assert(np < _np_quad[i_int]);
  v[0] = _table_quad[i_int][np].r;
  v[1] = _table_quad[i_int][np].s;
}


inline uint16 ElementIsoParamQUAD::nOfIntPoints() {
  return _np_quad[i_int];
}


inline double ElementIsoParamHEXAHEDRON::g_weight(uint16 np) {
  assert(np < _np_hexahedron[i_int]);
  assert(det.size() != 0);

  return _table_hexahedron[i_int][np].w * det[np];
}


inline void ElementIsoParamHEXAHEDRON::np2rst (uint16 np, double *v) {
  assert(np < _np_hexahedron[i_int]);
  v[0] = _table_hexahedron[i_int][np].r;
  v[1] = _table_hexahedron[i_int][np].s;
  v[2] = _table_hexahedron[i_int][np].t;
}


inline uint16 ElementIsoParamHEXAHEDRON::nOfIntPoints() {
  return _np_hexahedron[i_int];
}

} // namespace nla3d

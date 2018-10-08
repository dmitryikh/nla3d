// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"

namespace nla3d {

class ElementINTER0 : public ElementTWIN_VERTEX {
public:
  ElementINTER0 ();

  void pre();

  void buildK();

  void update();

  // stiffness
  double kn = 0.0, ks = 0.0;

  math::Vec<3> n; //local axis of spring

  math::Vec<3> strains; // displacement jump

  //postproc procedures
  bool getVector(math::Vec<3>* vector, vectorQuery code, uint16 gp, const double scale);
  bool getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale);
};

} //namespace nla3d

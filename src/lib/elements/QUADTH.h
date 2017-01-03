// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"
#include "elements/isoparametric.h"
#include "FEStorage.h"
#include "solidmech.h"

namespace nla3d {

//4-node 2D QUAD  element for steady thermal analysis
class ElementQUADTH : public ElementIsoParamQUAD {
  public:
    ElementQUADTH () {
      intOrder = 2;
      type = ElementType::QUADTH;
    }

    //solving procedures
    void pre();
    void buildK();
    void buildC();
    void buildM() { };
    void update();
    math::Mat<2,4> make_B (uint16 nPoint);  // make derivatives matrix

    // conductivity coef ( W/(K m), for example)
    double k = 0.0;
    double rho = 0.0; // density
    double c = 0.0; // capacity

    // volume flux
    double volFlux = 0.0;
};

class SurfaceLINETH : public ElementIsoParamLINE {
  public:
    SurfaceLINETH () {
      intOrder = 2;
      type = ElementType::SurfaceLINETH;
    }

    //solving procedures
    void pre();
    void buildK();
    void buildC() { };
    void buildM() { };
    void update();

    double flux = 0.0;
    double htc = 0.0;
    math::Vec<2> etemp = {0.0, 0.0};
};
} // nla3d namespace

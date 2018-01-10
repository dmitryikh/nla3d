// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/QUADTH.h"

namespace nla3d {
using namespace math;


void ElementQUADTH::pre() {
  if (det.size()==0) {
    makeJacob();
  }

  // register element equations
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::TEMP});
  }
}

void ElementQUADTH::buildK() {
  MatSym<4> Ke;  // element stiff. matrix
  Ke.zero();

  Vec<4> Fe; // rhs of element equations
  Fe.zero();

  MatSym<2> mat_k;
  mat_k.zero();
  mat_k.comp(0,0) = k;
  mat_k.comp(1,1) = k;

  // build Ke
  double dWt; //Gaussian quadrature weight
  for (uint16 np=0; np < nOfIntPoints(); np++) {
    dWt = intWeight(np);
    Mat<2,4> matB = make_B(np);
    matBTDBprod(matB, mat_k, dWt, Ke);

    // for volume flux
    if (volFlux != 0.0) {
      Fe += formFunc(np) * (volFlux * dWt);
    }
  }// loop over integration points

  assembleK(Ke, Fe, {Dof::TEMP});
}


void ElementQUADTH::buildC() {
  MatSym<4> Ce;  // element stiff. matrix
  Ce.zero();


  // build Ke
  double dWt; //Gaussian quadrature weight
  for (uint16 np=0; np < nOfIntPoints(); np++) {
    dWt = intWeight(np);
    Vec<4> ff = formFunc(np);
    for (uint16 i = 0; i < 4; i++)
      for (uint16 j = i; j < 4; j++) 
        Ce.comp(i, j) += ff[i] * ff[j] * rho * c * dWt;
  }// loop over integration points

  assembleC(Ce, {Dof::TEMP});
}

Mat<2,4> ElementQUADTH::make_B(uint16 np) {
  return Mat<2,4>(NiXj[np][0][0], NiXj[np][1][0], NiXj[np][2][0], NiXj[np][3][0],
                  NiXj[np][0][1], NiXj[np][1][1], NiXj[np][2][1], NiXj[np][3][1]);
}


void ElementQUADTH::update() {

}

void SurfaceLINETH::pre() {
  makeJacob();

  // register element equations
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::TEMP});
  }

  // if environment temperature isn't defined for node 2, we suppose that they are equal
  if (etemp[1] == 0.0) {
    etemp[1] = etemp[0];
  }
}


void SurfaceLINETH::buildK() {
  MatSym<2> Ke;  // element stiff. matrix
  Ke.zero();

  Vec<2> Fe; // rhs of element equations
  Fe.zero();

  double dWt; //Gaussian quadrature weight

  // flux over surface (line)
  if (flux != 0.0) {
    for (uint16 np = 0; np < nOfIntPoints(); np++) {
      dWt = intWeight(np);
      Vec<2> Ns = formFunc(np);
      Fe += Ns * (flux * dWt);
    }
  }

  // convection over surface (line)
  if (htc != 0.0) {
    for (uint16 np = 0; np < nOfIntPoints(); np++) {
      dWt = intWeight(np);
      Vec<2> Ns = formFunc(np);
      MatSym<2> tmp;
      tmp.comp(0,0) = Ns[0] * Ns[0] * htc * dWt;
      tmp.comp(0,1) = Ns[0] * Ns[1] * htc * dWt;
      tmp.comp(1,0) = tmp.comp(0,1);
      tmp.comp(1,1) = Ns[1] * Ns[1] * htc * dWt;
      //matAATprod(Ns, htc * dWt, tmp);
      Ke += tmp;
      matBVprod(tmp, etemp, 1.0, Fe);
    }
  }

  assembleK(Ke, Fe, {Dof::TEMP});
}

void SurfaceLINETH::update() {

}


} // namespace nla3d

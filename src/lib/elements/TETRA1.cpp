#include "elements/TETRA1.h"

namespace nla3d {

ElementTETRA1::ElementTETRA1 () {
  type = ElementType::TETRA1;
}

void ElementTETRA1::pre () {
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::TEMP});
  }
}

void ElementTETRA1::buildK() {
  math::Mat<4,4> matS(1. , storage->getNode(getNodeNumber(0)).pos[0] , storage->getNode(getNodeNumber(0)).pos[1] , storage->getNode(getNodeNumber(0)).pos[2] ,
                      1. , storage->getNode(getNodeNumber(1)).pos[0] , storage->getNode(getNodeNumber(1)).pos[1] , storage->getNode(getNodeNumber(1)).pos[2] ,
                      1. , storage->getNode(getNodeNumber(2)).pos[0] , storage->getNode(getNodeNumber(2)).pos[1] , storage->getNode(getNodeNumber(2)).pos[2] ,
                      1. , storage->getNode(getNodeNumber(3)).pos[0] , storage->getNode(getNodeNumber(3)).pos[1] , storage->getNode(getNodeNumber(3)).pos[2]);

  vol = matS.det()/6.;

  math::MatSym<4> matKe;
  matKe.zero();

  math::Mat<3,4> matB;
  matB.zero();

  math::MatSym<3> matC;
  matC.zero();

  makeC(matC);
  makeB(matB);
  math::matBTDBprod(matB, matC, vol, matKe);
  assembleK(matKe, {Dof::TEMP});
}

void ElementTETRA1::update () {
  math::Mat<3,4> matB;
  matB.zero();

  math::MatSym<3> matC;
  matC.zero();

  makeC(matC);
  makeB(matB);

  math::Vec<4> U;
  for (uint16 i = 0; i < getNNodes(); i++) {
    U[i] = storage->getNodeDofSolution(getNodeNumber(i), Dof::TEMP);
  }

  flux.zero();
  math::matBVprod(matB, U, k, flux);
}

void ElementTETRA1::makeB(math::Mat<3,4> &B)
{
  double *B_L = B.ptr();
  double b[4], c[4], d[4];

  int x=0, y = 1, z=2;

  double x12 = storage->getNode(getNodeNumber(0)).pos[x] - storage->getNode(getNodeNumber(1)).pos[x];
  double x13 = storage->getNode(getNodeNumber(0)).pos[x] - storage->getNode(getNodeNumber(2)).pos[x];
  double x14 = storage->getNode(getNodeNumber(0)).pos[x] - storage->getNode(getNodeNumber(3)).pos[x];
  double x23 = storage->getNode(getNodeNumber(1)).pos[x] - storage->getNode(getNodeNumber(2)).pos[x];
  double x24 = storage->getNode(getNodeNumber(1)).pos[x] - storage->getNode(getNodeNumber(3)).pos[x];
  double x34 = storage->getNode(getNodeNumber(2)).pos[x] - storage->getNode(getNodeNumber(3)).pos[x];

  double x21 = -1.*x12;
  double x31 = -1.*x13;
  double x32 = -1.*x23;
  double x42 = -1.*x24;
  double x43 = -1.*x34;

  double y12 = storage->getNode(getNodeNumber(0)).pos[y] - storage->getNode(getNodeNumber(1)).pos[y];
  double y13 = storage->getNode(getNodeNumber(0)).pos[y] - storage->getNode(getNodeNumber(2)).pos[y];
  double y14 = storage->getNode(getNodeNumber(0)).pos[y] - storage->getNode(getNodeNumber(3)).pos[y];
  double y23 = storage->getNode(getNodeNumber(1)).pos[y] - storage->getNode(getNodeNumber(2)).pos[y];
  double y24 = storage->getNode(getNodeNumber(1)).pos[y] - storage->getNode(getNodeNumber(3)).pos[y];
  double y34 = storage->getNode(getNodeNumber(2)).pos[y] - storage->getNode(getNodeNumber(3)).pos[y];

  double y21 = -1.*y12;
  double y31 = -1.*y13;
  double y32 = -1.*y23;
  double y42 = -1.*y24;
  double y43 = -1.*y34;

  double z12 = storage->getNode(getNodeNumber(0)).pos[z] - storage->getNode(getNodeNumber(1)).pos[z];
  double z13 = storage->getNode(getNodeNumber(0)).pos[z] - storage->getNode(getNodeNumber(2)).pos[z];
  double z14 = storage->getNode(getNodeNumber(0)).pos[z] - storage->getNode(getNodeNumber(3)).pos[z];
  double z23 = storage->getNode(getNodeNumber(1)).pos[z] - storage->getNode(getNodeNumber(2)).pos[z];
  double z24 = storage->getNode(getNodeNumber(1)).pos[z] - storage->getNode(getNodeNumber(3)).pos[z];
  double z34 = storage->getNode(getNodeNumber(2)).pos[z] - storage->getNode(getNodeNumber(3)).pos[z];

  double z21 = -1.*z12;
  double z31 = -1.*z13;
  double z32 = -1.*z23;
  double z42 = -1.*z24;
  double z43 = -1.*z34;

  b[0] = y42*z32 - y32*z42;
  b[1] = y31*z43 - y34*z13;
  b[2] = y24*z14 - y14*z24;
  b[3] = y13*z21 - y12*z31;

  c[0] = x32*z42-x42*z32;
  c[1] = x43*z31-x13*z34;
  c[2] = x14*z24-x24*z14;
  c[3] = x21*z13-x31*z12;

  d[0] = x42*y32-x32*y42;
  d[1] = x31*y43-x34*y13;
  d[2] = x24*y14-x14*y24;
  d[3] = x13*y21-x12*y31;

  const double A = -1./6./vol;
  B_L[0] = b[0]*A;
  B_L[1] = b[1]*A;
  B_L[2] = b[2]*A;
  B_L[3] = b[3]*A;
  B_L[4] = c[0]*A;
  B_L[5] = c[1]*A;
  B_L[6] = c[2]*A;
  B_L[7] = c[3]*A;
  B_L[8] = d[0]*A;
  B_L[9] = d[1]*A;
  B_L[10] = d[2]*A;
  B_L[11] = d[3]*A;
}

void ElementTETRA1::makeC (math::MatSym<3> &C) {
  C.comp(0,0) = -k;
  C.comp(1,1) = -k;
  C.comp(2,2) = -k;
}

bool ElementTETRA1::getScalar(double* scalar, scalarQuery query, uint16 gp, const double scale) {
  if (query == scalarQuery::VOL){
     *scalar += vol*scale;
    return true;
  }
  return false;
}

bool  ElementTETRA1::getVector(math::Vec<3>* vector, vectorQuery query, uint16 gp, const double scale) {
  switch (query) {
    case vectorQuery::FLUX:
      *vector += flux*scale;
      return true;
    case vectorQuery::GRADT:
      *vector += flux*(scale/k);
      return true;
  }  
  return false;
}

} //namespace nla3d

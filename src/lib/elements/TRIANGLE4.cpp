#include "elements/TRIANGLE4.h"

#include <iostream>

using namespace std;

namespace nla3d {

ElementTRIANGLE4::ElementTRIANGLE4 () {
  type = ElementType::TRIANGLE4;
  state = PlaneState::Stress;
}

void ElementTRIANGLE4::pre () {
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY});
  }
}

// here stiffness matrix is built
void ElementTRIANGLE4::buildK() {
  // Ke will store element stiffness matrix in global coordinates
  math::MatSym<6> matKe;
  matKe.zero();

  // matB is strain matrix
  math::Mat<3,6> matB;
  matB.zero();

  // matC is 3d elastic  matrix
  math::MatSym<3> matC;
  matC.zero();
  //only for area 
  Eigen::MatrixXd matS(3,3);
  matS.setZero();
  matS << 1. , storage->getNode(getNodeNumber(0)).pos[0] , storage->getNode(getNodeNumber(0)).pos[1] ,
          1. , storage->getNode(getNodeNumber(1)).pos[0] , storage->getNode(getNodeNumber(1)).pos[1] ,
          1. , storage->getNode(getNodeNumber(2)).pos[0] , storage->getNode(getNodeNumber(2)).pos[1];
  area = matS.determinant()/2.;
  // fill here matC
  makeC(matC);
  // fill here matB
  makeB(matB);
 
  math::matBTDBprod(matB, matC, area, matKe);
  // start assemble procedure. Here we should provide element stiffness matrix and an order of 
  // nodal DoFs in the matrix.
  assembleK(matKe, {Dof::UX, Dof::UY});
}

// after solution it's handy to calculate stresses, strains and other stuff in elements.
void ElementTRIANGLE4::update () {
  // matB is strain matrix
  math::Mat<3,6> matB;
  matB.zero();

  // matC is 3d elastic  matrix
  math::MatSym<3> matC;
  matC.zero();
  // fill here matC
  makeC(matC);
  // fill here matB
  makeB(matB);
  // get nodal solutions from storage
  math::Vec<6> U;
  for (uint16 i = 0; i < getNNodes(); i++) {
    U[i*2 + 0] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U[i*2 + 1] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
  }
  // restore strains
  math::matBVprod(matB, U, 1.0, strains);
  // restore stresses
  math::matBVprod(matC, strains, 1.0, stress);
}

void ElementTRIANGLE4::makeB(math::Mat<3,6> &B)
{
    double *B_L = B.ptr();
    double b[3], c[3];
    b[0] = storage->getNode(getNodeNumber(1)).pos[1] - storage->getNode(getNodeNumber(2)).pos[1];
    b[1] = storage->getNode(getNodeNumber(2)).pos[1] - storage->getNode(getNodeNumber(0)).pos[1];
    b[2] = storage->getNode(getNodeNumber(0)).pos[1] - storage->getNode(getNodeNumber(1)).pos[1];
    c[0] = storage->getNode(getNodeNumber(2)).pos[0] - storage->getNode(getNodeNumber(1)).pos[0];
    c[1] = storage->getNode(getNodeNumber(0)).pos[0] - storage->getNode(getNodeNumber(2)).pos[0];
    c[2] = storage->getNode(getNodeNumber(1)).pos[0] - storage->getNode(getNodeNumber(0)).pos[0];

    const double A = 1./2./area;

    B_L[0*6 + 0] = b[0]*A;
    B_L[0*6 + 2] = b[1]*A;
    B_L[0*6 + 4] = b[2]*A;
    B_L[1*6 + 1] = c[0]*A;
    B_L[1*6 + 3] = c[1]*A;
    B_L[1*6 + 5] = c[2]*A;
    B_L[2*6 + 0] = c[0]*A;
    B_L[2*6 + 1] = b[0]*A;
    B_L[2*6 + 2] = c[1]*A;
    B_L[2*6 + 3] = b[1]*A;
    B_L[2*6 + 4] = c[2]*A;
    B_L[2*6 + 5] = b[2]*A;
}

void ElementTRIANGLE4::makeC (math::MatSym<3> &C) {
    //Plane deformed state
    if (state == PlaneState::Strain){
      const double A = E*(1.-my)/((1.+my)*(1.-2.*my));
      C.comp(0,0) = 1.*A;
      C.comp(0,1) = my/(1-my)*A;
      C.comp(1,1) = 1.*A;
      C.comp(2,2) = (1.-2.*my)/(2.*(1.-my))*A;
    }
    //Plane stress state
    if (state == PlaneState::Stress){
      const double A = E*(1.-my*my);
      C.comp(0,0) = 1.*A;
      C.comp(0,1) = my*A; 
      C.comp(1,1) = 1.*A;
      C.comp(2,2) = (1.-my)/2.*A; 
    }
}

} //namespace nla3d

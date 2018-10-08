#include "elements/TETRA0.h"

using namespace std;

namespace nla3d {

ElementTETRA0::ElementTETRA0 () {
  type = ElementType::TETRA0;
}

void ElementTETRA0::pre () {
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY, Dof::UZ});
  }
}

// here stiffness matrix is built
void ElementTETRA0::buildK() {
  Eigen::MatrixXd matS(4,4);
  matS.setZero();
  matS<< 1. , storage->getNode(getNodeNumber(0)).pos[0] , storage->getNode(getNodeNumber(0)).pos[1] , storage->getNode(getNodeNumber(0)).pos[2] ,
          1. , storage->getNode(getNodeNumber(1)).pos[0] , storage->getNode(getNodeNumber(1)).pos[1] , storage->getNode(getNodeNumber(1)).pos[2] ,
          1. , storage->getNode(getNodeNumber(2)).pos[0] , storage->getNode(getNodeNumber(2)).pos[1] , storage->getNode(getNodeNumber(2)).pos[2] ,
          1. , storage->getNode(getNodeNumber(3)).pos[0] , storage->getNode(getNodeNumber(3)).pos[1] , storage->getNode(getNodeNumber(3)).pos[2];

  vol = matS.determinant()/6.;
  // Ke will store element stiffness matrix in global coordinates
  math::MatSym<12> matKe;
  matKe.zero();

  // matB is strain matrix
  math::Mat<6,12> matB;
  matB.zero();

  // matC is 3d elastic  matrix
  math::MatSym<6> matC;
  matC.zero();

  // fill here matC
  makeC(matC);
  // fill here matB
  makeB(matB);  

  math::matBTDBprod(matB, matC, vol, matKe);

  if ((alpha != 0. && T != 0.) || strains.qlength() != 0. || stress.qlength() != 0.){
    //node forces calculations
    math::Vec<12> Fe;
    Fe.zero();

    math::Mat<12,6> matBTC;
    matBTC = matB.transpose()*matC.toMat();

    //mechanical initial stress
    if (stress.qlength() != 0.){
      strains.zero();
      math::Mat<6,6> matP;
      matP = matC.toMat().inv(matC.toMat().det());
      strains = matP*stress;
    }
    
    //termal initial strains
    if (alpha != 0. && T != 0.){
      //temp node forces
      math::Vec<6> tStrains = {alpha*T,alpha*T,alpha*T,0.,0.,0.};
      strains = strains + tStrains;
    }


    //mechanical initial strains
    math::matBVprod(matBTC, strains, -vol, Fe);

    assembleK(matKe, Fe, {Dof::UX, Dof::UY, Dof::UZ});
  }
  else{
    assembleK(matKe, {Dof::UX, Dof::UY, Dof::UZ});
  }
}

// after solution it's handy to calculate stresses, strains and other stuff in elements.
void ElementTETRA0::update () {
  // matB is strain matrix
  math::Mat<6,12> matB;
  matB.zero();

  // matC is 3d elastic  matrix
  math::MatSym<6> matC;
  matC.zero();

  // fill here matC
  makeC(matC);
  // fill here matB
  makeB(matB);
  // get nodal solutions from storage
  math::Vec<12> U;
  for (uint16 i = 0; i < getNNodes(); i++) {
    U[i*3 + 0] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U[i*3 + 1] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
    U[i*3 + 2] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UZ);
  }
  
  //restore strains
  strains.zero();
  math::matBVprod(matB, U, -1.0, strains);
  
  stress.zero();
  math::matBVprod(matC, strains, 1.0, stress);
}

void ElementTETRA0::makeB(math::Mat<6,12> &B)
{
    double *B_L = B.ptr();
    double b[4], c[4], d[4];
    //Eigen::MatrixXd mb(3,3), mc(3,3), md(3,3);
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
    for (int i = 0; i < 4; i++){
      B_L[0 + 3*i] = b[i]*A;
      B_L[13 + 3*i] = c[i]*A;
      B_L[26 + 3*i] = d[i]*A;
      B_L[36 + 3*i] = c[i]*A;
      B_L[37 + 3*i] = b[i]*A;
      B_L[49 + 3*i] = d[i]*A;
      B_L[50 + 3*i] = c[i]*A;
      B_L[60 + 3*i] = d[i]*A;
      B_L[62 + 3*i] = b[i]*A;
    }
}

void ElementTETRA0::makeC (math::MatSym<6> &C) {
  const double A = E/((1.+my)*(1.-2.*my));

  C.comp(0,0) = (1.-my)*A;
  C.comp(0,1) = my*A;
  C.comp(0,2) = my*A;
  C.comp(1,1) = (1.-my)*A;
  C.comp(1,2) = my*A;
  C.comp(2,2) = (1.-my)*A;

  C.comp(3,3) = (1./2.-my)*A;
  C.comp(4,4) = (1./2.-my)*A;
  C.comp(5,5) = (1./2.-my)*A;
}

bool ElementTETRA0::getScalar(double* scalar, scalarQuery query, uint16 gp, const double scale) {
  if (query == scalarQuery::VOL){
     *scalar += vol;
    return true;
  }
  return false;
}

bool  ElementTETRA0::getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale) {
  if (query == tensorQuery::C){
      tensor->comp(0,0) += strains[0];
      tensor->comp(1,1) += strains[1];
      tensor->comp(2,2) += strains[2];
      tensor->comp(0,1) += strains[3];
      tensor->comp(1,2) += strains[4];
      tensor->comp(0,2) += strains[5];
      return true;
  }
  if (query == tensorQuery::E){
    tensor->comp(0,0) += stress[0];
    tensor->comp(1,1) += stress[1];
    tensor->comp(2,2) += stress[2];
    tensor->comp(0,1) += stress[3];
    tensor->comp(1,2) += stress[4];
    tensor->comp(0,2) += stress[5];
    return true;
  }
  if (query == tensorQuery::TSTRAIN){
    tensor->comp(0,0) += alpha*T;
    tensor->comp(1,1) += alpha*T;
    tensor->comp(2,2) += alpha*T;
    tensor->comp(0,1) += 0.;
    tensor->comp(1,2) += 0.;
    tensor->comp(0,2) += 0.;
    return true;
  }
  
  return false;
}
} //namespace nla3d

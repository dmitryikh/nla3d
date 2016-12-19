#include "elements/TETRA.h"

namespace nla3d {

ElementTETRA::ElementTETRA () {
  Element::number_of_nodes = 4;
  number_of_dimensions = 3;
  nodes = new uint32[Element::number_of_nodes];

}

void ElementTETRA::pre () {
  for (uint16 i = 0; i < Element::n_nodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), Dof::UX);
    storage->addNodeDof(getNodeNumber(i), Dof::UY);
    storage->addNodeDof(getNodeNumber(i), Dof::UZ);
  }
}

// here stiffness matrix is built
void ElementTETRA::build () {
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

  double vol = 1.0; // TODO: calculate REAL!! volume here

  math::matBTDBprod(matB, matC, vol, matKe);

  // start assemble procedure. Here we should provide element stiffness matrix and an order of 
  // nodal DoFs in the matrix.
  assemble(matKe, {Dof::UX, Dof::UY, Dof::UZ});
}

// after solution it's handy to calculate stresses, strains and other stuff in elements.
void ElementTETRA::update () {
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
  for (uint16 i = 0; i < Element::n_nodes(); i++) {
    U[i*3 + 0] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U[i*3 + 1] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
    U[i*3 + 2] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UZ);
  }

  // restore strains
  math::matBVprod(matB, U, 1.0, strains);
  
  // restore stresses
  // convert from symmetric mat to full mat, because func matBVprod now works only with full mat
  // здесь не работает.. Это связано с тем, что сейчас два класса матриц mat и mat2.. Можно
  // использовать eigen и забить на свои классы, но eigen на функциях matBVprod и  matBTDBprod
  // работают медленней собственных кодов.. поэтому приходится для критических мест держать свои
  // матрицы и вектора. Хорошо бы иметь трансляторы из одного класса в другой
  // math::matBVprod(matC.toMat(), strains, 1.0, strains);
}

void ElementTETRA::makeB(math::Mat<6,12> &B)
{
    double *B_L = B.ptr();
    for (uint16 i=0; i < 4; i++) {
        // need to write correct values here!!!
        B_L[0*12+(i*3+0)] += 1;//exx
        B_L[1*12+(i*3+0)] += 1;//exy
        B_L[1*12+(i*3+1)] += 1;//exy
        B_L[2*12+(i*3+0)] += 1;//exz
        B_L[2*12+(i*3+2)] += 1;//exz
        B_L[3*12+(i*3+1)] += 1;//eyy
        B_L[4*12+(i*3+1)] += 1;//eyz
        B_L[4*12+(i*3+2)] += 1;//eyz
        B_L[5*12+(i*3+2)] += 1;//ezz
    }
}

void ElementTETRA::makeC (math::MatSym<6> &C) {
    // fill matrix here..
}

} //namespace nla3d

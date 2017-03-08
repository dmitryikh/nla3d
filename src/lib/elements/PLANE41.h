// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"
#include "elements/isoparametric.h"
#include "FEStorage.h"
#include "solidmech.h"

namespace nla3d {

//4-node 2D QUAD nonlinear element based on mixed approach
class ElementPLANE41 : public ElementIsoParamQUAD {
  public:
    ElementPLANE41 () {
      intOrder = 2;
      type = ElementType::PLANE41;
    }

    ElementPLANE41 (const ElementPLANE41& from) {
      operator=(from);
    }

    //solving procedures
    void pre();
    void buildK();
    void update();
    math::Mat<3,8> make_B (uint16 nPoint);  //функция создает линейную матрицу [B]
    math::Mat<4,8> make_Bomega (uint16 nPoint); //функция создает линейную матрицу [Bomega]

    //postproc procedures
    bool getScalar(double* scalar, scalarQuery code, uint16 gp, const double scale);
    bool getVector(math::Vec<3>* vector, vectorQuery code, uint16 gp, const double scale);
    bool getTensor(math::MatSym<3>* tensor, tensorQuery code, uint16 gp, const double scale);

    // internal element data
    // S[0] - Sx  S[1] - Sy S[2] - Sxy
    std::vector<math::Vec<3> > S; //S[номер т. интегр.][номер напряжения] - напряжения Пиолы-Кирхгоффа
    // C[0] - C11 C[1] - C22  C[2] - C12
    std::vector<math::Vec<3> > C; //C[номер т. интегр.][номер деформ.] - компоненты матрицы меры деформации
    // O[0] - dU/dx O[1] - dU/dy  O[2] - dV/dx  O[3] - dV/dy
    std::vector<math::Vec<4> > O; //S[номер т. интегр.][номер омеги]

    // addition data
    static const solidmech::tensorComponents components[3];
    static const uint16 num_components;


    template <uint16 el_dofs_num>
    void assembleK(const math::Mat<el_dofs_num,el_dofs_num> &Ke, const math::Vec<el_dofs_num> &Qe);
};

template <uint16 el_dofs_num>
void ElementPLANE41::assembleK(const math::Mat<el_dofs_num,el_dofs_num> &Ke, const math::Vec<el_dofs_num> &Qe)
{
  uint16 dim = 2;
  uint16 eds = 8;
  uint16 eldofs = 1;
  Dof::dofType nodeDofVec[] = {Dof::UX, Dof::UY};
  Dof::dofType elementDofVec[] = {Dof::HYDRO_PRESSURE};
  assert(getNNodes() * dim + eldofs == el_dofs_num);
  // upper diagonal process for nodal dofs
  for (uint16 i=0; i < getNNodes(); i++)
    for (uint16 j=i; j < getNNodes(); j++)
      for (uint16 di=0; di < dim; di++)
        for (uint16 dj=0; dj < dim; dj++)
          if ((i==j) && (dj>di)) continue;
          else
            storage->addValueK(nodes[i], nodeDofVec[di],nodes[j],nodeDofVec[dj], Ke[i*dim+di][j*dim+dj]);

  //upper diagonal process for nodes-el dofs
  for (uint16 i=0; i < getNNodes(); i++)
    for(uint16 di=0; di < dim; di++) {
      uint32 noEq = storage->getNodeDofEqNumber(nodes[i], nodeDofVec[di]);
      for (uint16 dj=0; dj < eldofs; dj++) {
        uint32 elEq = storage->getElementDofEqNumber(getElNum(), elementDofVec[dj]);
        storage->addValueK(noEq, elEq, Ke[i*dim+di][eds+dj]);
      }
    }
  //upper diagonal process for el-el dofs
  for (uint16 di=0; di < eldofs; di++)
    for (uint16 dj=di; dj < eldofs; dj++) {
      uint32 elEqi = storage->getElementDofEqNumber(getElNum(), elementDofVec[di]);
      uint32 elEqj = storage->getElementDofEqNumber(getElNum(), elementDofVec[dj]);
      storage->addValueK(elEqi, elEqj, Ke[eds+di][eds+dj]);
    }

  for (uint16 i=0; i < getNNodes(); i++)
    for (uint16 di=0; di < dim; di++)
      storage->addValueF(nodes[i], nodeDofVec[di], Qe[i*dim+di]);

  for (uint16 di=0; di < eldofs; di++) {
    uint32 elEq = storage->getElementDofEqNumber(getElNum(), elementDofVec[di]);
    storage->addValueF(elEq, Qe[eds+di]);
  }
}

} // namespace nla3d 

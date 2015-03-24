// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"
#include "elements/element_lagrange.h"
#include "FEStorage.h"
#include "solidmech.h"

namespace nla3d {

//-------------------------------------------------------
//-------------------ElementPLANE41----------------------
//-------------------------------------------------------
//4-node 2D QUAD nonlinear element based on mixed approach
class ElementPLANE41 : public Element_PLANE4, public Element_Lagrange_Formulation<2,4> {
public:
	ElementPLANE41 () {
    Node::registerDofType(Dof::UX);
    Node::registerDofType(Dof::UY);
    Element::registerDofType(Dof::HYDRO_PRESSURE);
	}
	ElementPLANE41 (const ElementPLANE41& from) {
		operator=(from);
	}

  //solving procedures
	void pre();
	void build();
	void update ();
	math::Mat<3,8> make_B (uint16 nPoint);	//функция создает линейную матрицу [B]
	math::Mat<4,8> make_Bomega (uint16 nPoint);	//функция создает линейную матрицу [Bomega]

  //postproc procedures
	void getScalar(double& scalar, query::scalarQuery code, uint16 gp, const double scale);
	void getVector(double* vector, query::vectorQuery code, uint16 gp, const double scale);
	void getTensor(math::MatSym<3>& tensor, query::tensorQuery code, uint16 gp, const double scale);

  // internal element data
	// S[0] - Sx	S[1] - Sy	S[2] - Sxy
	std::vector<math::Vec<3> > S; //S[номер т. интегр.][номер напряжения] - напряжения Пиолы-Кирхгоффа
	// C[0] - C11	C[1] - C22	C[2] - C12
	std::vector<math::Vec<3> > C; //C[номер т. интегр.][номер деформ.] - компоненты матрицы меры деформации
	// O[0] - dU/dx	O[1] - dU/dy	O[2] - dV/dx	O[3] - dV/dy
	std::vector<math::Vec<4> > O; //S[номер т. интегр.][номер омеги]

  // addition data
	static const solidmech::tensorComponents components[3];
	static const uint16 num_components;


	template <uint16 el_dofs_num>
	void assemble (const math::Mat<el_dofs_num,el_dofs_num> &Ke, const math::Vec<el_dofs_num> &Qe);
};

template <uint16 el_dofs_num>
void ElementPLANE41::assemble (const math::Mat<el_dofs_num,el_dofs_num> &Ke, const math::Vec<el_dofs_num> &Qe)
{
  uint16 dim = 2;
  uint16 eds = 8;
  Dof::dofType nodeDofVec[] = {Dof::UX, Dof::UY};
  Dof::dofType elementDofVec[] = {Dof::HYDRO_PRESSURE};
	assert(Element::n_nodes() * dim + Element::n_dofs() == el_dofs_num);
	// upper diagonal process for nodal dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 j=i; j < Element::n_nodes(); j++)
			for (uint16 di=0; di < dim; di++)
				for (uint16 dj=0; dj < dim; dj++)
					if ((i==j) && (dj>di)) continue;
					else
						storage->Kij_add(nodes[i], nodeDofVec[di],nodes[j],nodeDofVec[dj], Ke[i*dim+di][j*dim+dj]);
	//upper diagonal process for nodes-el dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for(uint16 di=0; di < dim; di++)
			for (uint16 dj=0; dj < Element::n_dofs(); dj++)
				storage->Kij_add(nodes[i], nodeDofVec[di], -(int32)getElNum(), elementDofVec[dj], Ke[i*dim+di][eds+dj]);
	//upper diagonal process for el-el dofs
	for (uint16 di=0; di < Element::n_dofs(); di++)
		for (uint16 dj=di; dj < Element::n_dofs(); dj++)
			storage->Kij_add(-(int32)getElNum(), elementDofVec[di], -(int32)getElNum(), elementDofVec[dj],  Ke[eds+di][eds+dj]);

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < dim; di++)
			storage->Fi_add(nodes[i], nodeDofVec[di], Qe[i*dim+di]);

	for (uint16 di=0; di < Element::n_dofs(); di++)
		storage->Fi_add(-(int32)getElNum(), elementDofVec[di], Qe[eds+di]);
}

} // namespace nla3d 

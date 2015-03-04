#pragma once
#include "elements/element.h"
#include "elements/element_lagrange.h"
#include "FE_Storage.h"

//-------------------------------------------------------
//-------------------ElementPLANE41----------------------
//-------------------------------------------------------
//4-node 2D QUAD nonlinear element based on mixed approach
class ElementPLANE41 : public Element_PLANE4, public Element_Lagrange_Formulation<2,4> {
public:
	ElementPLANE41 ()
	{
		change_node_dofs_num(2, UX, UY);
		change_el_dofs_num(1,HYDRO_PRESSURE);
	}
	ElementPLANE41 (const ElementPLANE41& from)
	{
		operator=(from);
	}

  //solving procedures
	void pre();
	void build();
	void update ();
	Mat<3,8> make_B (uint16 nPoint);	//функция создает линейную матрицу [B]
	Mat<4,8> make_Bomega (uint16 nPoint);	//функция создает линейную матрицу [Bomega]

  //postproc procedures
	void getScalar(double& scalar, el_component code, uint16 gp, const double scale);
	void getTensor(MatSym<3>& tensor, el_tensor code, uint16 gp, const double scale);

  // internal element data
	// S[0] - Sx	S[1] - Sy	S[2] - Sxy
	vector<Vec<3> > S; //S[номер т. интегр.][номер напряжения] - напряжения Пиолы-Кирхгоффа
	// C[0] - C11	C[1] - C22	C[2] - C12
	vector<Vec<3> > C; //C[номер т. интегр.][номер деформ.] - компоненты матрицы меры деформации
	// O[0] - dU/dx	O[1] - dU/dy	O[2] - dV/dx	O[3] - dV/dy
	vector<Vec<4> > O; //S[номер т. интегр.][номер омеги]

  // addition data
	static const tensorComponents components[3];
	static const uint16 num_components;


	template <uint16 el_dofs_num>
	void assemble (const Mat<el_dofs_num,el_dofs_num> &Ke, const Vec<el_dofs_num> &Qe);
};

template <uint16 el_dofs_num>
void ElementPLANE41::assemble (const Mat<el_dofs_num,el_dofs_num> &Ke, const Vec<el_dofs_num> &Qe)
{
	uint16 eds = Element::n_nodes()*Node::n_dofs(); // el's dofs start number
	assert(el_dofs_num == eds+Element::n_dofs());
	// upper diagonal process for nodal dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 j=i; j < Element::n_nodes(); j++)
			for (uint16 di=0; di < Node::n_dofs(); di++)
				for (uint16 dj=0; dj < Node::n_dofs(); dj++)
					if ((i==j) && (dj>di)) continue;
					else
						storage->Kij_add(nodes[i],di,nodes[j],dj, Ke[i*Node::n_dofs()+di][j*Node::n_dofs()+dj]);
	//upper diagonal process for nodes-el dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for(uint16 di=0; di < Node::n_dofs(); di++)
			for (uint16 dj=0; dj < Element::n_dofs(); dj++)
				storage->Kij_add(nodes[i],di, -(int32)getElNum(), dj, Ke[i*Node::n_dofs()+di][eds+dj]);
	//upper diagonal process for el-el dofs
	for (uint16 di=0; di < Element::n_dofs(); di++)
		for (uint16 dj=di; dj < Element::n_dofs(); dj++)
			storage->Kij_add(-(int32)getElNum(), di, -(int32)getElNum(), dj,  Ke[eds+di][eds+dj]);

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
			storage->Fi_add(nodes[i],di, Qe[i*Node::n_dofs()+di]);

	for (uint16 di=0; di < Element::n_dofs(); di++)
		storage->Fi_add(-(int32)getElNum(), di, Qe[eds+di]);
}


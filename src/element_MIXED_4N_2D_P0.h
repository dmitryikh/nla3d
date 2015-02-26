#pragma once

#include "Element.h"
#include "Element_Lagrange_Formulation.h"

//-------------------------------------------------------
//-------------------MIXED_4N_2D_P0----------------------
//-------------------------------------------------------
//4-node 2D QUAD nonlinear element based on mixed approach
class MIXED_4N_2D_P0 : public Element_4N_2D, public Element_Lagrange_Formulation<2,4>
{
public:
	MIXED_4N_2D_P0 ()
	{
		change_node_dofs_num(2, UX, UY);
		change_el_dofs_num(1,HYDRO_PRESSURE);
	}
	MIXED_4N_2D_P0 (const MIXED_4N_2D_P0& from)
	{
		operator=(from);
	}

  //solving procedures
	void pre (uint32 el, FE_Storage_Interface *storage);
	void build (uint32 el, FE_Storage_Interface *storage);
	void update (uint32 el, FE_Storage_Interface *storage);
	Mat<3,8> make_B (uint16 nPoint);	//функция создает линейную матрицу [B]
	Mat<4,8> make_Bomega (uint16 nPoint);	//функция создает линейную матрицу [Bomega]

  //postproc procedures
	double getComponent (uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage);
  void getTensor (uint16 gp, el_tensor code, uint32 el, FE_Storage_Interface *storage, MatSym<3> &tensor);

  // internal element data
	// S[0] - Sx	S[1] - Sy	S[2] - Sxy
	vector<Vec<3>> S; //S[номер т. интегр.][номер напряжения] - напряжения Пиолы-Кирхгоффа
	// C[0] - C11	C[1] - C22	C[2] - C12
	vector<Vec<3>> C; //C[номер т. интегр.][номер деформ.] - компоненты матрицы меры деформации
	// O[0] - dU/dx	O[1] - dU/dy	O[2] - dV/dx	O[3] - dV/dy
	vector<Vec<4>> O; //S[номер т. интегр.][номер омеги]

  // addition data
	static const mat_comp components[3];
	static const uint16 num_components;
};


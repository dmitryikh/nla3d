#pragma once

#include "Element.h"
#include "Element_Lagrange_Formulation.h"


//4-node 2D QUAD linear element (pure displacement approach)
class DISP_4N_2D : public Element_4N_2D, public Element_Lagrange_Formulation<2,4>
{
public:
	DISP_4N_2D ()
	{
		change_node_dofs_num(2, UX, UY);
		change_el_dofs_num(0);
	}
	DISP_4N_2D (const DISP_4N_2D& from)
	{
		operator=(from);
	}
	void pre (uint32 el, FE_Storage_Interface *storage);
	void build (uint32 el, FE_Storage_Interface *storage);
	void update (uint32 el, FE_Storage_Interface *storage);
	virtual double getComponent (uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage);
	Mat<3,8> make_B (uint16 nPoint);	//функция создает линейную матрицу [B]
};

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
	void pre (uint32 el, FE_Storage_Interface *storage);
	void build (uint32 el, FE_Storage_Interface *storage);
	void update (uint32 el, FE_Storage_Interface *storage);
	virtual double getComponent (uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage);

	Mat<3,8> make_B (uint16 nPoint);	//функция создает линейную матрицу [B]
	Mat<4,8> make_Bomega (uint16 nPoint);	//функция создает линейную матрицу [Bomega]

	vector<Vec<3>> S; //S[номер т. интегр.][номер напряжения] - напряжения Пиолы-Кирхгоффа
	// S[0] - Sx	S[1] - Sy	S[2] - Sxy
	vector<Vec<3>> C; //C[номер т. интегр.][номер деформ.] - компоненты матрицы меры деформации
	// C[0] - C11	C[1] - C22	C[2] - C12
	vector<Vec<4>> O; //S[номер т. интегр.][номер омеги]
	// O[0] - dU/dx	O[1] - dU/dy	O[2] - dV/dx	O[3] - dV/dy
};

//8-node brick linear element based on mixed approach
class MIXED_8N_3D_P0 : public Element_8N_3D, public Element_Lagrange_Formulation<3,8>
{
public:
	MIXED_8N_3D_P0 ()
	{
		change_node_dofs_num(3, UX, UY, UZ);
		change_el_dofs_num(1,HYDRO_PRESSURE);
	}
	MIXED_8N_3D_P0 (const MIXED_8N_3D_P0& from)
	{
		operator=(from);
	}
	void pre (uint32 el, FE_Storage_Interface *storage);
	void build (uint32 el, FE_Storage_Interface *storage);
	void update (uint32 el, FE_Storage_Interface *storage);
	virtual double getComponent (uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage);
	MIXED_8N_3D_P0& operator= (const MIXED_8N_3D_P0& from)
	{
		Element_Lagrange_Formulation<3,8>::operator= (from);
		Element_8N_3D::operator= (from);
		return *this;
	}
private:
	Mat<6,24> make_b_d (uint16 nPoint);
	Mat<1,24> make_b_v (uint16 nPoint);
};


//8-node brick nonlinear element based on mixed approach
class MIXED_8N_3D_P0_NL : public Element_8N_3D, public Element_Lagrange_Formulation<3,8>
{
public:
	MIXED_8N_3D_P0_NL ()
	{
		change_node_dofs_num(3, UX, UY, UZ);
		change_el_dofs_num(1,HYDRO_PRESSURE);
	}
	MIXED_8N_3D_P0_NL (const MIXED_8N_3D_P0_NL& from)
	{
		operator=(from);
	}
	void pre (uint32 el, FE_Storage_Interface *storage);
	void build (uint32 el, FE_Storage_Interface *storage);
	void update (uint32 el, FE_Storage_Interface *storage);
	virtual double getComponent (uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage);

	Mat<6,24> make_B_L (uint16 nPoint);	//функция создает линейную матрицу [B]
	Mat<9,24> make_B_NL (uint16 nPoint);	//функция создает линейную матрицу [Bomega]
	Mat<9,9> make_S (uint16 nPoint);
	Mat<6,9> make_Omega (uint16 nPoint);
	vector<Vec<6>> S; //S[номер т. интегр.][номер напряжения] - напряжения Пиолы-Кирхгоффа
	// S[0] - Sx	S[1] - Sy	S[2] - Sz	S[3] - Sxy	S[4] - Syz	S[5] - Sxz
	vector<Vec<6>> C; //C[номер т. интегр.][номер деформ.] - компоненты тензора меры деформации
	// C[0] - C11	C[1] - C22	C[2] - C33	C[3] - C12	C[4] - C23	C[5] - C[5]
	vector<Vec<9>> O; //S[номер т. интегр.][номер омеги]
	// O[0]-dU/dx	O[1]-dU/dy	O[2]-dU/dz	O[3]-dV/dx	O[4]-dV/dy	O[5]-dV/dz	O[6]-dW/dx	O[7]-dW/dy	O[8]-dW/dz
};

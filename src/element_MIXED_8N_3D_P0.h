#pragma once

#include "Element.h"
#include "Element_Lagrange_Formulation.h"

//-------------------------------------------------------
//-------------------MIXED_8N_3D_P0----------------------
//-------------------------------------------------------
//8-node brick nonlinear element based on mixed approach
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

  //solving procedures
	void pre ();
	void build ();
	void update ();

	void make_B_L (uint16 nPoint, Mat2<6,24> &B);	//функция создает линейную матрицу [B]
	void make_B_NL (uint16 nPoint,  Mat2<9,24> &B); //функция создает линейную матрицу [Bomega]
	void make_S (uint16 nPoint, MatSym<9> &B);
	void make_Omega (uint16 nPoint, Mat2<6,9> &B);

  //postproc procedures
	void getScalar(double& scalar, el_component code, uint16 gp, const double scale);
	void getTensor(MatSym<3>& tensor, el_tensor code, uint16 gp, const double scale);

  // internal element data
	// S[0] - Sx	S[1] - Sy	S[2] - Sz	S[3] - Sxy	S[4] - Syz	S[5] - Sxz
	vector<Vec<6>> S; //S[номер т. интегр.][номер напряжения] - напряжения Пиолы-Кирхгоффа
	// C[0] - C11	C[1] - C22	C[2] - C33	C[3] - C12	C[4] - C23	C[5] - C[5]
	vector<Vec<6>> C; //C[номер т. интегр.][номер деформ.] - компоненты тензора меры деформации
	// O[0]-dU/dx	O[1]-dU/dy	O[2]-dU/dz	O[3]-dV/dx	O[4]-dV/dy	O[5]-dV/dz	O[6]-dW/dx	O[7]-dW/dy	O[8]-dW/dz
	vector<Vec<9>> O; //S[номер т. интегр.][номер омеги]

  // addition data
	static const mat_comp components[6];
	static const uint16 num_components;
};

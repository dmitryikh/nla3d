#pragma once
#include "elements/element.h"
#include "elements/element_lagrange.h"
#include "FE_Storage.h"
#include <Eigen/Dense>

//-------------------------------------------------------
//-------------------ElementSOLID81----------------------
//-------------------------------------------------------
//8-node brick nonlinear element based on mixed approach
class ElementSOLID81 : public Element_SOLID8, public Element_Lagrange_Formulation<3,8> {
public:
	ElementSOLID81 () {
		change_node_dofs_num(3, UX, UY, UZ);
		change_el_dofs_num(1,HYDRO_PRESSURE);
	}
	ElementSOLID81 (const ElementSOLID81& from)
	{
		operator=(from);
	}

  //solving procedures
	void pre();
	void build();
	void update();

	void make_B_L (uint16 nPoint, Eigen::Ref<Eigen::MatrixXd> B);	//функция создает линейную матрицу [B]
	void make_B_NL (uint16 nPoint, Eigen::Ref<Eigen::MatrixXd> B); //функция создает линейную матрицу [Bomega]
	void make_S (uint16 nPoint, Eigen::Ref<Eigen::MatrixXd> SMat);
	void make_Omega (uint16 nPoint, Eigen::Ref<Eigen::MatrixXd> Omega);

  //postproc procedures
	void getScalar(double& scalar, el_component code, uint16 gp, const double scale);
    void getVector(double* vector, el_vector code, uint16 gp, const double scale);
	void getTensor(MatSym<3>& tensor, el_tensor code, uint16 gp, const double scale);

  // internal element data
	//S[M_XX], S[M_XY], S[M_XZ], S[M_YY], S[M_YZ], S[M_ZZ]
	vector<Vec<6> > S; //S[номер т. интегр.][номер напряжения] - напряжения Пиолы-Кирхгоффа
	//C[M_XX], C[M_XY], C[M_XZ], C[M_YY], C[M_YZ], C[M_ZZ]
	vector<Vec<6> > C; //C[номер т. интегр.][номер деформ.] - компоненты тензора меры деформации
	// O[0]-dU/dx	O[1]-dU/dy	O[2]-dU/dz	O[3]-dV/dx	O[4]-dV/dy	O[5]-dV/dz	O[6]-dW/dx	O[7]-dW/dy	O[8]-dW/dz
	vector<Vec<9> > O; //S[номер т. интегр.][номер омеги]

  // addition data
	static const tensorComponents components[6];
	static const uint16 num_components;


	template <uint16 dimM, uint16 dimN>
	void assemble2(MatSym<dimM> &Kuu, Mat2<dimM,dimM> &Kup, Mat2<dimN,dimN> &Kpp, Vec<dimM> &Fu, Vec<dimN> &Fp);
	template <uint16 dimM>
	void assemble3(MatSym<dimM> &Kuu, Vec<dimM> &Kup, double Kpp, Vec<dimM> &Fu, double Fp);
  void assemble4(Eigen::Ref<Eigen::MatrixXd> Ke, Eigen::Ref<Eigen::MatrixXd> Fe);
};

template <uint16 dimM, uint16 dimN>
void ElementSOLID81::assemble2(MatSym<dimM> &Kuu, Mat2<dimM,dimM> &Kup, Mat2<dimN,dimN> &Kpp, Vec<dimM> &Fu, Vec<dimN> &Fp) 
{
	assert (Element::n_nodes()*Node::n_dofs() == dimM);
	assert (Element::n_dofs() == dimN);
	double *Kuu_p = Kuu.ptr();
	double *Kup_p = Kup.ptr();
	double *Kpp_p = Kpp.ptr();
	double *Fu_p = Fu.ptr();
	double *Fp_p = Fp.ptr();                               

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
			for (uint16 j=i; j < Element::n_nodes(); j++)
				for (uint16 dj=di; dj < Node::n_dofs(); dj++)
				{
						storage->Kij_add(nodes[i],di,nodes[j],dj, *Kuu_p);
						Kuu_p++;
				}
	//upper diagonal process for nodes-el dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for(uint16 di=0; di < Node::n_dofs(); di++)
			for (uint16 dj=0; dj < Element::n_dofs(); dj++)
			{
				storage->Kij_add(nodes[i],di, -(int32)getElNum(), dj, *Kup_p);
				Kup_p++;
			}
	//upper diagonal process for el-el dofs
	for (uint16 di=0; di < Element::n_dofs(); di++)
		for (uint16 dj=di; dj < Element::n_dofs(); dj++)
		{
			storage->Kij_add(-(int32)getElNum(), di, -(int32)getElNum(), dj,  *Kpp_p);
			Kpp_p++;
		}

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
		{
			storage->Fi_add(nodes[i],di, *Fu_p);
			Fu_p++;
		}

	for (uint16 di=0; di < Element::n_dofs(); di++)
	{
		storage->Fi_add(-(int32)getElNum(), di, *Fp_p);
		Fp_p++;
	}
}

template <uint16 dimM>
void ElementSOLID81::assemble3(MatSym<dimM> &Kuu, Vec<dimM> &Kup, double Kpp, Vec<dimM> &Fu, double Fp) 
{
	assert (Element::n_nodes()*Node::n_dofs() == dimM);
	assert (Element::n_dofs() == 1);
	double *Kuu_p = Kuu.ptr();
	double *Kup_p = Kup.ptr();
	double *Fu_p = Fu.ptr();

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
		for (uint16 j=i; j < Element::n_nodes(); j++)
			
				for (uint16 dj=0; dj < Node::n_dofs(); dj++)
				{
					if ((i==j) && (dj<di)) continue;
					else
					{
						storage->Kij_add(nodes[i],di,nodes[j],dj, *Kuu_p);
						Kuu_p++;
					}
				}
	//upper diagonal process for nodes-el dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for(uint16 di=0; di < Node::n_dofs(); di++)
			for (uint16 dj=0; dj < Element::n_dofs(); dj++)
			{
				storage->Kij_add(nodes[i],di, -(int32)getElNum(), dj, *Kup_p);
				Kup_p++;
			}
	//upper diagonal process for el-el dofs
			storage->Kij_add(-(int32)getElNum(), 0, -(int32)getElNum(), 0,  Kpp);

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
		{
			storage->Fi_add(nodes[i],di, *Fu_p);
			Fu_p++;
		}
		storage->Fi_add(-(int32)getElNum(), 0, Fp);
}


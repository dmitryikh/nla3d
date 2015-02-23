#pragma once
#include "sys.h"
#include "Element.h"

// dim - смысл размерности геометрии элемента
// квадратичные фф??? надо реализовать
template <uint32 dim, uint32 nodes_num>
class Element_Lagrange_Formulation
{
public:
	static vector<Vec<nodes_num>> gN; // calculated form function on ecah int point
	vector<Mat<dim,dim>> Jacob;	//Jacob inv matrix
	vector<Mat<dim,nodes_num>> NjXi; //derivates form function / local coordinates
	vector<double>	det;	//Jacobian
	const static Mat<4,4> gPoint; //integration points & weights for Gauss integration
	const static Mat<4,4> gWeight; // supported schemes 1x1x1 to 4x4x4
	void make_Jacob(uint32 el, Node** nodes);	//function to calculate all this staff above
	void npoint_to_ii(uint16 nPoint, uint16 *ii); //by number of gauss point find indexes
	void npoint_to_xi(uint16 nPoint, double *xi); //by number of gauss point find local coordinates
	double g_weight(uint16 nPoint);
	const static int16 _signs[8][3];
};


template <uint32 dim, uint32 nodes_num> const Mat<4,4> Element_Lagrange_Formulation<dim,nodes_num>::gPoint = Mat<4,4>(0.0f, 0.0f, 0.0f, 0.0f,
												-0.5773502691896f, +0.5773502691896f, 0.0f, 0.0f,
												-0.7745966692415f, 0.0f, +0.7745966692415f, 0.0f,
												-0.8611363115941f, -0.3399810435849f, +0.3399810435849f, +0.8611363115941f);
template <uint32 dim, uint32 nodes_num> const Mat<4,4> Element_Lagrange_Formulation<dim,nodes_num>::gWeight = Mat<4,4>(2.0f, 0.0f, 0.0f, 0.0f,
												1.0f, 1.0f, 0.0f, 0.0f,
												0.5555555555556f, 0.8888888888889f, 0.5555555555556f, 0.0f,
												0.3478548451375f, 0.6521451548625f, 0.6521451548625f, 0.3478548451375f);
template <uint32 dim, uint32 nodes_num> vector<Vec<nodes_num>> Element_Lagrange_Formulation<dim,nodes_num>::gN;

template <uint32 dim, uint32 nodes_num>const  int16 Element_Lagrange_Formulation<dim,nodes_num>::_signs[8][3] = {-1,-1,-1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,1,1,-1,1,1,1,1,-1,1,1};


template <uint32 dim, uint32 nodes_num>
void Element_Lagrange_Formulation<dim,nodes_num>::npoint_to_xi(uint16 nPoint, double *xi)
{
	uint16 ii[dim];
	npoint_to_ii(nPoint,ii);
	xi[0] = 0.0;
	xi[1] = 0.0;
	xi[2] = 0.0;
	for (uint16 i = 0; i < dim; i++)
			xi[i]=gPoint[Element::n_int()-1][ii[i]];
}

template <uint32 dim, uint32 nodes_num>
void Element_Lagrange_Formulation<dim,nodes_num>::npoint_to_ii(uint16 nPoint, uint16 *ii)
{
	uint16 residual=0;
	for (uint16 i=0; i<dim; i++)
	{
		ii[dim-1-i]= (uint16) (nPoint-residual)/npow(Element::n_int(),dim-i-1);
		residual += npow(Element::n_int(),dim-i-1)*ii[dim-1-i];
	}
}

template <uint32 dim, uint32 nodes_num> 
void Element_Lagrange_Formulation<dim,nodes_num>::make_Jacob(uint32 el, Node** nodes) //is el need?
{
	double xi[3]; // local coordinates (xi, eta, zeta)
	uint16 nPoint = 0; // current int point
	//once fill up matrix of calculated function form 
	if (gN.size()==0)
	{
		gN.assign(npow(Element::n_int(), dim), Vec<nodes_num>());
		for (nPoint=0; (int32) nPoint < npow(Element::n_int(),dim); nPoint++)
		{
			npoint_to_xi(nPoint, xi);
			for (uint16 i = 0; i < nodes_num; i++)
				gN[nPoint][i] = 1.0/nodes_num*(1+_signs[i][0]*xi[0])*(1+_signs[i][1]*xi[1])*(1+_signs[i][2]*xi[2]);
		}
	}

	Jacob.clear();
	Jacob.assign(npow(Element::n_int(),dim), Mat<dim,dim>());
	det.clear();
	det.assign(npow(Element::n_int(),dim), 0.0);
	NjXi.clear();
	NjXi.assign(npow(Element::n_int(),dim), Mat<dim,nodes_num>());
	double inv_det;
	Mat<dim,nodes_num> dN; //производные функции формы
	Mat<dim,dim> J;
	for (nPoint=0; (int32) nPoint < npow(Element::n_int(),dim); nPoint++)
	{
		for (uint16 der = 0; der < dim; der++)
		{
			npoint_to_xi(nPoint, xi);
			xi[der] = 0.0;
			for (uint16 nod = 0; nod < nodes_num; nod++)
				dN[der][nod] = _signs[nod][der]*1.0/nodes_num*(1+_signs[nod][0]*xi[0])*(1+_signs[nod][1]*xi[1])*(1+_signs[nod][2]*xi[2]);
		}
		J.zero();
		for (uint16 i = 0; i < dim; i++)
			for (uint16 j = 0; j < dim; j++)
				for (uint16 nod = 0; nod < nodes_num; nod++)	
					J[i][j] += dN[i][nod]*nodes[nod]->pos[j];

		det[nPoint] = J.det(); // determinant of Jacob matrix // 
		if (det[nPoint] < 1.0e-20)  // check for geometry form error 
			warning("Element_form::make_Jacob: det too small in element %d", el); 
		// обращение матрицы Якоби
		inv_det = 1.0/det[nPoint];
		Jacob[nPoint] = J.inv(det[nPoint]);
		// производные функций формы по глоб. координатам
		for (uint16 nod = 0; nod < nodes_num; nod++)
			for (uint16 i=0; i < dim; i++)
				for (uint16 j=0; j< dim; j++)
					NjXi[nPoint][i][nod]+=Jacob[nPoint][i][j]*dN[j][nod];
	}

}
template <uint32 dim, uint32 nodes_num> 
double Element_Lagrange_Formulation<dim,nodes_num>::g_weight(uint16 nPoint)
{
	assert(det.size() != 0);
	uint16 ii[3];
	npoint_to_ii(nPoint, ii);
	double res = det[nPoint];
	for (uint16 i=0; i < dim; i++)
		res *= gWeight[Element::n_int()-1][ii[i]];
	return res;
}
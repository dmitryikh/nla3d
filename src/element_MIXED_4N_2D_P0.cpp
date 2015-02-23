#include "element_MIXED_4N_2D_P0.h"

const mat_comp MIXED_4N_2D_P0::components[3] = {M_XX, M_YY, M_XY};
const uint16 MIXED_4N_2D_P0::num_components = 3;

//------------------MIXED_4N_2D_P0--------------------
void MIXED_4N_2D_P0::pre (uint32 el, FE_Storage_Interface *storage)
{
	S.assign(npow(n_int(),n_dim()), Vec<3>(0.0f, 0.0f, 0.0f));
	C.assign(npow(n_int(),n_dim()), Vec<3>(1.0f, 1.0f, 0.0f));
	O.assign(npow(n_int(),n_dim()), Vec<4>(0.0f, 0.0f, 0.0f, 0.0f));
	if (det.size()==0) 
	{
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->element_nodes(el, nodes_p);	
		make_Jacob(el, nodes_p);
	}
}
void MIXED_4N_2D_P0::build (uint32 el, FE_Storage_Interface *storage)
{
//построение матрицы жесткости и вектора нагрузки для элемента
	Mat<8,8> Kuu; //матрица жесткости элемента
	Mat<8,1> Kup;
	Mat<9,9> Ke;
	double Kpp = 0.0;
	double Fp = 0.0;
	Vec<8> Qe; //вектор узловых сил элемента
	Vec<9> Fe;
//M_XX =  0,
//M_XY =  1,
//M_XZ =  2,
//M_YY =  3,
//M_YZ =  4,
//M_ZZ =  5
  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
	MatSym<3> matD_d;
	Vec<3> vecD_p;
	double p_e = storage->get_qi_n(-(int32)el, 0);
	GlobStates.setdouble(States::DOUBLE_HYDPRES, p_e);
  Material *mat = (Material*) GlobStates.getptr(States::PTR_CURMATER);
	double k = mat->getK0();
	double dWt; //множитель при суммировании квадратур Гаусса
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		dWt = g_weight(nPoint);
    //TODO: we need to reimplement this functions, because now Materaial class is different
    // all meterial functions are waiting [C] for 3D case. So we need to use CVec here.
    CVec[M_XX] = C[nPoint][0];
    CVec[M_YY] = C[nPoint][1];
    CVec[M_XY] = C[nPoint][2];
		mat->getDdDp_UP(num_components, components, CVec.ptr(), matD_d.ptr(), vecD_p.ptr());
    //matD_d will be 3x3 symmetric matrix. We need to convert it onto 3x3 usual matrix
    Mat<3,3> matE_c = matD_d.toMat();
    Mat<3,1> matE_p;
    matE_p[0][0] = vecD_p[0];
    matE_p[1][0] = vecD_p[1];
    matE_p[2][0] = vecD_p[2];
		double J = Material::getJ(CVec.ptr());

		//Mat<3,3> matE_c = storage->getMaterial().mat_E_c(ANALYSIS_2D_PLANE_STRAIN, C[nPoint].ptr(), p_e).toMat<3,3>();
		//Mat<3,1> matE_p = storage->getMaterial().mat_E_p(ANALYSIS_2D_PLANE_STRAIN, C[nPoint].ptr()).toMat<3,1>();
		Mat<3,8> matB = make_B(nPoint);
		//матрица S для матричного умножения
		Mat<4,4> matS = Mat<4,4>(S[nPoint][0],S[nPoint][2], 0.0, 0.0, 
								S[nPoint][2],S[nPoint][1], 0.0, 0.0,
								0.0, 0.0, S[nPoint][0],S[nPoint][2],
								0.0, 0.0, S[nPoint][2],S[nPoint][1]);
		//матрица Омега.используется для составления 
		//матр. накопленных линейных деформаций к текущему шагу
		Mat<3,4> matO = Mat<3,4>(O[nPoint][0], 0.0, O[nPoint][2], 0.0,
								0.0, O[nPoint][1], 0.0, O[nPoint][3],
								O[nPoint][1], O[nPoint][0], O[nPoint][3], O[nPoint][2]);

		Mat<4,8> matBomega = make_Bomega(nPoint);
		Mat<3,8> matBl = matO * matBomega;
		matB += matBl;
    //TODO: somehow here is 2.0 multiplayer.. I need to figure out why..
		Kuu += (matB.transpose() * matE_c * matB * 2.0 + matBomega.transpose() * matS * matBomega)*dWt;
		//double J = (1.0+O[nPoint][0])*(1.0+O[nPoint][3])-O[nPoint][1]*O[nPoint][2];
		Fp += (J - 1 - p_e/k)*dWt;
		Qe += (matB.transpose() * S[nPoint] * dWt);
		Kup+= matB.transpose()*matE_p *dWt;
		Kpp -= 1.0/k*dWt;

	}//прошлись по всем точкам интегрирования
	
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);

	//сборка в одну матрицу
	for (uint16 i=0; i < 8; i++)
		for (uint16 j=0; j < 8; j++)
			Ke[i][j] = Kuu[i][j];
	for (uint16 i=0; i<8; i++)
		Ke[i][8] = Kup[i][0];
	for (uint16 i=0; i<8; i++)
		Ke[8][i] = Kup[i][0];
	Ke[8][8] = Kpp;
	for (uint16 i=0; i < 8; i++)
		Fe[i] = -Qe[i];
	Fe[8] = -Fp;
	//загнать в глоб. матрицу жесткости и узловых сил
	assemble(el, Ke, Fe, storage);
}
//
inline Mat<3,8> MIXED_4N_2D_P0::make_B (uint16 nPoint) {
	Mat<3,8> B = Mat<3,8>(NjXi[nPoint][0][0], 0.0f, NjXi[nPoint][0][1], 0.0f, NjXi[nPoint][0][2], 0.0f, NjXi[nPoint][0][3], 0.0f,
									0.0f, NjXi[nPoint][1][0], 0.0f, NjXi[nPoint][1][1], 0.0f, NjXi[nPoint][1][2], 0.0f, NjXi[nPoint][1][3],
									NjXi[nPoint][1][0], NjXi[nPoint][0][0], NjXi[nPoint][1][1], NjXi[nPoint][0][1], NjXi[nPoint][1][2], NjXi[nPoint][0][2], NjXi[nPoint][1][3], NjXi[nPoint][0][3]);
	return B;
}
//
Mat<4,8> MIXED_4N_2D_P0::make_Bomega (uint16 nPoint)
{
	Mat<4,8> Bomega = Mat<4,8>(NjXi[nPoint][0][0], 0.0f, NjXi[nPoint][0][1], 0.0f, NjXi[nPoint][0][2], 0.0f, NjXi[nPoint][0][3], 0.0f,
										NjXi[nPoint][1][0], 0.0f, NjXi[nPoint][1][1], 0.0f, NjXi[nPoint][1][2], 0.0f, NjXi[nPoint][1][3], 0.0f,
										0.0f, NjXi[nPoint][0][0], 0.0f, NjXi[nPoint][0][1], 0.0f, NjXi[nPoint][0][2], 0.0f, NjXi[nPoint][0][3],
										0.0f, NjXi[nPoint][1][0], 0.0f, NjXi[nPoint][1][1], 0.0f, NjXi[nPoint][1][2], 0.0f, NjXi[nPoint][1][3]);
	return Bomega;
}
//
void MIXED_4N_2D_P0::update (uint32 el, FE_Storage_Interface *storage)
{
	Vec<9> U; //вектор перемещений элемента
	// получаем вектор перемещений элемента из общего решения
	storage->get_q_e(el, U.ptr());
	Vec<8> Un;
  Material *mat = (Material*) GlobStates.getptr(States::PTR_CURMATER);
  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
	for (uint16 i=0;i<8;i++)
		Un[i] = U[i];
	//восстанавливаем преращение давления
	double p_e = U[8];;
	GlobStates.setdouble(States::DOUBLE_HYDPRES, p_e);

  for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++) {
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		Mat<4,8> matBomega = make_Bomega(nPoint);
		O[nPoint] = matBomega * Un;
		C[nPoint][0] = 1.0f + 2*O[nPoint][0]+1.0f*(O[nPoint][0]*O[nPoint][0]+O[nPoint][2]*O[nPoint][2]);
		C[nPoint][1]=1.0f + 2*O[nPoint][3]+1.0f*(O[nPoint][3]*O[nPoint][3]+O[nPoint][1]*O[nPoint][1]);
		C[nPoint][2]=O[nPoint][1]+O[nPoint][2]+O[nPoint][0]*O[nPoint][1]+O[nPoint][2]*O[nPoint][3];
		//восстановление напряжений Пиолы-Кирхгоффа из текущего состояния
    //TODO: Now Material class is differ
		//S[nPoint] = storage->getMaterial().getS(ANALYSIS_2D_PLANE_STRAIN, C[nPoint].ptr(), p_e).toVec<3>();
    // all meterial functions are waiting [C] for 3D case. So we need to use CVec here.
    CVec[M_XX] = C[nPoint][0];
    CVec[M_YY] = C[nPoint][1];
    CVec[M_XY] = C[nPoint][2];
		mat->getS_UP (num_components, components, CVec.ptr(), S[nPoint].ptr());
	}
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);
}

double MIXED_4N_2D_P0::getComponent(uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage)
{
	//////////////////////////////////////
	//see codes in sys.h
	//gp - needed gauss point 

	assert (gp < n_int()*n_int());
	Mat<3,3> matX;
	Mat<3,3> matS;
	Mat<3,3> matT;
	Vec<3> E_vec;
	double res, J;
	matX.zero();
	matX[0][0] = 1+O[gp][0];
	matX[0][1] = O[gp][1];
	matX[1][0] = O[gp][2];
	matX[1][1] = 1+O[gp][3];
	matX[2][2] = 1;

	E_vec[0] = 1.0/2.0*(C[gp][0]-1);
	E_vec[1] = 1.0/2.0*(C[gp][1]-1);
	E_vec[2] = 1.0f/2.0*C[gp][2];

	matS[0][0] = S[gp][0];
	matS[0][1] = S[gp][2];
	matS[1][0] = S[gp][2];
	matS[1][1] = S[gp][1];
  //TODO: actually we have Szz..

	J= (matX[0][0]*matX[1][1] - matX[0][1]*matX[1][0])*matX[2][2];
	double J2 = sqrt(C[gp][0]*C[gp][1]-C[gp][2]*C[gp][2]);
	matT = matX * matS * matX.transpose() * (1.0f/J);


	switch (code)
	{
	case E_X:
		res = E_vec[0];
		break;
	case E_Y:
		res = E_vec[1];
		break;
	case E_XY:
		res = E_vec[2];
		break;
	case E_Z:
		res = 0; //нет такой компоненты в плоской задаче ПДС
		break;
	case E_VOL:
		res = J2;
		break;
	case S_X:
		res = matT[0][0];
		break;
	case S_Y:
		res = matT[1][1];
		break;
	case S_XY:
		res = matT[0][1];
		break;
	case S_Z:
		res = matT[2][2];
		break;
	case S_P:
		res =  storage->get_qi_n(-(int32)el, 0);
		break;
	default:
    //TODO: for 2D and 3D elements here are different components available..
    //some how we need to write into vtk file only relevant components
		//warning("Element::getComponent: error in code %d", code);
		res=0.0;
	}
	return res;
}
// TODO: implement this
Mat<3,3> MIXED_4N_2D_P0::getTensor (uint16 gp, el_tensor code, uint32 el, FE_Storage_Interface *storage) {
	Vec<3> cmass(0.0, 0.0, 0.0);
	Vec<3> xi;
	for (uint16 i = 0; i < Element::n_nodes(); i ++)
	{
		storage->get_node_pos(storage->getElement(el).node_num(i), xi.ptr(), true);
		cmass = cmass + xi;
	}
	cmass = cmass * (1.0/Element::n_nodes());
	cmass[1]=0.0;
	//TODO: Storage::el_center(el, def/undef)
	cmass = cmass * (1.0/cmass.lenght());
	
	Vec<3> nr(cmass[0],cmass[1],0.0);
	Vec<3> nz(0.0,0.0,1.0);
	Vec<3> nt(cmass[1],-cmass[0],0.0);
	
	Mat<3,3> Q(cmass[0],cmass[2],0.0,
				0.0,0.0,1.0,
				cmass[2],	-cmass[0],	0.0);

	Mat<3,3> matX;
	Mat<3,3> matS;
	Mat<3,3> matT;

	matX.zero();
	matX[0][0] = 1+O[gp][0];
	matX[0][1] = O[gp][1];
	matX[1][0] = O[gp][2];
	matX[1][1] = 1+O[gp][3];
	matX[2][2] = 1;

	matS[0][0] = S[gp][0];
	matS[0][1] = S[gp][2];
	matS[0][2] = 0.0;
	matS[1][0] = S[gp][2];
	matS[1][1] = S[gp][1];
	matS[1][2] = 0.0;
	matS[2][0] = 0.0;
	matS[2][1] = 0.0;
	matS[2][2] = 0.0;//it's worng. Actually we have Szz
	matT = matX * matS * matX.transpose() * (1.0f/matX.det());




	switch (code)
	{
	case TENS_COUCHY:
		return matT;
		break;
  default:
    error ("What kind of Tensor you want?");
	}

	return Mat<3,3>(0.0,0.0,0.0,
					0.0,0.0,0.0,
					0.0,0.0,0.0);
}

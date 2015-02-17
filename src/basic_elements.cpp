#include "basic_elements.h"

//------------------DISP_4N_2D--------------------
void DISP_4N_2D::pre (uint32 el, FE_Storage_Interface *storage)
{
	if (det.size()==0) 
	{
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->element_nodes(el, nodes_p);	
		make_Jacob(el, nodes_p);
	}
}
void DISP_4N_2D::build (uint32 el, FE_Storage_Interface *storage)
{
	//построение матрицы жесткости и вектора нагрузки для элемента
	Mat<8,8> Ke; //матрица жесткости элемента
	Vec<8> Fe;
	double G = storage->getMaterial().getGxy0();
	double mu = storage->getMaterial().getMuxy0();
	double lambda = 2*mu*G/(1-2*mu);
	Mat<3,3> D(2*G+lambda,lambda,0.0,
				lambda, 2*G+lambda,0.0,
				0.0,	0.0,	G);
	Mat<3,8> B;

	double dWt; //множитель при суммировании квадратур Гаусса
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		dWt = g_weight(nPoint);
		B = make_B(nPoint);
		Ke += (B.transpose() * D * B)*dWt;	//TODO: реализовать B^T*D*B - как отдельный метод Mat<i,j> для быстрого перемножения
	}//прошлись по всем точкам интегрирования
	assemble(el, Ke, Fe, storage);
}

inline Mat<3,8> DISP_4N_2D::make_B (uint16 nPoint) {
	Mat<3,8> B = Mat<3,8>(NjXi[nPoint][0][0], 0.0f, NjXi[nPoint][0][1], 0.0f, NjXi[nPoint][0][2], 0.0f, NjXi[nPoint][0][3], 0.0f,
									0.0f, NjXi[nPoint][1][0], 0.0f, NjXi[nPoint][1][1], 0.0f, NjXi[nPoint][1][2], 0.0f, NjXi[nPoint][1][3],
									NjXi[nPoint][1][0], NjXi[nPoint][0][0], NjXi[nPoint][1][1], NjXi[nPoint][0][1], NjXi[nPoint][1][2], NjXi[nPoint][0][2], NjXi[nPoint][1][3], NjXi[nPoint][0][3]);
	return B;
}

void DISP_4N_2D::update (uint32 el, FE_Storage_Interface *storage)
{
}

double DISP_4N_2D::getComponent(uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage)
{
	//see codes in sys.h
	//gp - needed gauss point 
	//TODO: написать функцию!
	return 0.0;
	//assert (gp < n_int()*n_int());
	//Mat<3,3> matX;
	//Mat<3,3> matS;
	//Mat<3,3> matT;
	//Vec<3> E_vec;
	//double res, J;
	//matX.zero();
	//matX[0][0] = 1+O[gp][0];
	//matX[0][1] = O[gp][1];
	//matX[1][0] = O[gp][2];
	//matX[1][1] = 1+O[gp][3];
	//matX[2][2] = 1;

	//E_vec[0] = 1.0/2.0*(C[gp][0]-1);
	//E_vec[1] = 1.0/2.0*(C[gp][1]-1);
	//E_vec[2] = 1.0f/2.0*C[gp][2];

	//matS[0][0] = S[gp][0];
	//matS[0][1] = S[gp][2];
	//matS[1][0] = S[gp][2];
	//matS[1][1] = S[gp][1];

	//J= (matX[0][0]*matX[1][1] - matX[0][1]*matX[1][0])*matX[2][2];
	//double J2 = sqrt(C[gp][0]*C[gp][1]-C[gp][2]*C[gp][2]);
	//matT = matX * matS * matX.transpose() * (1.0f/J);

	////double pos_x, pos_y;
	////pos_x = 0.0; pos_y = 0.0;
	////Vec<3> npos;
	////for (uint16 i = 0; i < Element::nNodes(); i++)
	////{
	////	storage->get_node_pos(nodes[i], npos.ptr());
	////	pos_x += gN[gp][i]*npos[0];
	////	pos_y += gN[gp][i]*npos[1];
	////}


	//switch (code)
	//{
	//case E_X:
	//	res = E_vec[0];
	//	break;
	//case E_Y:
	//	res = E_vec[1];
	//	break;
	//case E_XY:
	//	res = E_vec[2];
	//	break;
	//case E_Z:
	//	res = 0; //нет такой компоненты в плоской задаче ПДС
	//	break;
	//case E_VOL:
	//	res = J2;
	//	break;
	//case S_X:
	//	res = matT[0][0];
	//	break;
	//case S_Y:
	//	res = matT[1][1];
	//	break;
	//case S_XY:
	//	res = matT[0][1];
	//	break;
	//case S_Z:
	//	res = matT[2][2];
	//	break;
	//case S_P:
	//	res =  storage->get_qi_n(-(int32)el, 0);
	//	break;
	////case POS_X:
	////	res = pos_x;
	////	break;
	////case POS_Y:
	////	res = pos_y;
	////	break;
	//default:
	//	//warning("Element::getComponent: error in code %d", code);
	//	res=0.0;
	//}
	//return res;
}

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
	double p_e = storage->get_qi_n(-(int32)el, 0);

	double dWt; //множитель при суммировании квадратур Гаусса
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		dWt = g_weight(nPoint);
		Mat<3,3> matE_c = storage->getMaterial().mat_E_c(ANALYSIS_2D_PLANE_STRAIN, C[nPoint].ptr(), p_e).toMat<3,3>();
		Mat<3,1> matE_p = storage->getMaterial().mat_E_p(ANALYSIS_2D_PLANE_STRAIN, C[nPoint].ptr()).toMat<3,1>();
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
		Kuu += (matB.transpose() * matE_c * matB *2 + matBomega.transpose() * matS * matBomega)*dWt;
		double J = (1.0+O[nPoint][0])*(1.0+O[nPoint][3])-O[nPoint][1]*O[nPoint][2];
		Fp += (J - 1 - p_e/storage->getMaterial().getK0())*dWt;
		Qe += (matB.transpose() * S[nPoint] * dWt);
		Kup+= matB.transpose()*matE_p*dWt;
		Kpp -= 1.0/storage->getMaterial().getK0()*dWt;

	}//прошлись по всем точкам интегрирования
	
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
	for (uint16 i=0;i<8;i++)
		Un[i] = U[i];
	//восстанавливаем преращение давления
	double p_e = U[8];;

		for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		Mat<4,8> matBomega = make_Bomega(nPoint);
		O[nPoint] = matBomega * Un;
		C[nPoint][0] = 1.0f + 2*O[nPoint][0]+1.0f*(O[nPoint][0]*O[nPoint][0]+O[nPoint][2]*O[nPoint][2]);
		C[nPoint][1]=1.0f + 2*O[nPoint][3]+1.0f*(O[nPoint][3]*O[nPoint][3]+O[nPoint][1]*O[nPoint][1]);
		C[nPoint][2]=O[nPoint][1]+O[nPoint][2]+O[nPoint][0]*O[nPoint][1]+O[nPoint][2]*O[nPoint][3];
		//восстановление напряжений Пиолы-Кирхгоффа из текущего состояния
		S[nPoint] = storage->getMaterial().getS(ANALYSIS_2D_PLANE_STRAIN, C[nPoint].ptr(), p_e).toVec<3>();
	}
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
		//warning("Element::getComponent: error in code %d", code);
		res=0.0;
	}
	return res;
}


//-----------------MIXED_8N_3D_P0-----------------
void MIXED_8N_3D_P0::pre (uint32 el, FE_Storage_Interface *storage)
{
	if (det.size()==0) 
	{
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->element_nodes(el, nodes_p);
		make_Jacob(el, nodes_p);
		delete[] nodes_p;
	}
}

void MIXED_8N_3D_P0::build (uint32 el, FE_Storage_Interface *storage)
{
	//построение матрицы жесткости и вектора нагрузки для элемента
	Mat<24,24> Kuu; //
	Mat<24,1> Kup;
	double Kpp = 0.0;
	Mat<25,25> Ke; //

	Vec<24> Qu; //вектор узловых сил элемента
	Vec<25> Fe;
	
	double dWt; //множитель при суммировании квадратур Гаусса
	Mat<6,6> D_d = storage->getMaterial().mat_D_d(ANALYSIS_3D).toMat<6,6>();
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		dWt = g_weight(nPoint);
		Mat<6,24> B_d = make_b_d(nPoint);
		Mat<1,24> B_v = make_b_v(nPoint);  
		Kuu += (B_d.transpose() * D_d * B_d)*dWt;
		Kup += B_v.transpose()*dWt;
		Kpp += 1.0/storage->getMaterial().getK0()*dWt;
	}

	//сборка в одну матрицу 
	//TODO: сделать эту процедуру более эффективной 
	//(например сразу работать с матрицей 25x25 и записывать результат только в нужные столбцы)
	for (uint16 i=0; i < 24; i++)
		for (uint16 j=0; j < 24; j++)
			Ke[i][j] = Kuu[i][j];
	for (uint16 i=0; i<24; i++)
		Ke[i][24] = Kup[i][0];
	for (uint16 i=0; i<24; i++)
		Ke[24][i] = Kup[i][0];
	Ke[24][24] = -Kpp;
	for (uint16 i=0; i < 24; i++)
		Fe[i] = Qu[i];

	assemble(el, Ke, Fe, storage);
}

inline Mat<6,24> MIXED_8N_3D_P0::make_b_d (uint16 nPoint)
{
	Mat<6,24> b_d;
	for (uint16 j = 0; j < 8; j++)
	{
		b_d[0][j*3+0]=2.0/3.0*NjXi[nPoint][0][j];  b_d[0][j*3+1]=-1.0/3.0*NjXi[nPoint][1][j];  b_d[0][j*3+2]=-1.0/3.0*NjXi[nPoint][2][j];

		b_d[1][j*3+0]=-1.0/3.0*NjXi[nPoint][0][j];  b_d[1][j*3+1]=2.0/3.0*NjXi[nPoint][1][j];  b_d[1][j*3+2]=-1.0/3.0*NjXi[nPoint][2][j];

		b_d[2][j*3+0]=-1.0/3.0*NjXi[nPoint][0][j];  b_d[2][j*3+1]=-1.0/3.0*NjXi[nPoint][1][j];  b_d[2][j*3+2]=2.0/3.0*NjXi[nPoint][2][j];

		b_d[3][j*3+0]=NjXi[nPoint][1][j];           b_d[3][j*3+1]=NjXi[nPoint][0][j];           b_d[3][j*3+2]=0.0;

		b_d[4][j*3+0]=0.0;                          b_d[4][j*3+1]=NjXi[nPoint][2][j];           b_d[4][j*3+2]=NjXi[nPoint][1][j];

		b_d[5][j*3+0]=NjXi[nPoint][2][j];           b_d[5][j*3+1]=0.0;                          b_d[5][j*3+2]=NjXi[nPoint][0][j];
	}
	return b_d;
}

inline Mat<1,24> MIXED_8N_3D_P0::make_b_v (uint16 nPoint)
{
	Mat<1,24> b_v;
	for (uint16 j = 0; j < 8; j++)
	{
		b_v[0][j*3+0]=NjXi[nPoint][0][j];  b_v[0][j*3+1]=NjXi[nPoint][1][j];  b_v[0][j*3+2]=NjXi[nPoint][2][j];
	}
	return b_v;
}

void MIXED_8N_3D_P0::update (uint32 el, FE_Storage_Interface *storage)
{
}

double MIXED_8N_3D_P0::getComponent(uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage)
{
	//see codes in sys.h
	//gp - needed gauss point 
	assert ((int32) gp < npow(Element::n_int(),3));
	Vec<6> vecEps;
	Vec<6> vecSig;
	Vec<6> vecEps_d;
	Vec<6> vecKroneker(1.0,1.0,1.0,0.0,0.0,0.0);
	Mat<6,6> D_d = storage->getMaterial().mat_D_d(ANALYSIS_3D).toMat<6,6>();
	double res, J;

	Vec<25> q_e; //вектор решения для степеней свободы узлов эл. и элемента
	storage->get_q_e(el, q_e.ptr()); // получаем вектор перемещений элемента из общего решения
	Vec<24> u_e;
	for (uint16 i=0; i < 24; i++)
		u_e[i] = q_e[i];
	vecEps_d = make_b_d(gp)*u_e;
	vecEps = vecEps_d + vecKroneker*(q_e[24]/(3*storage->getMaterial().getK0()));
	vecSig = D_d*vecEps_d+vecKroneker*q_e[24];
	J = vecEps[0]+vecEps[1]+vecEps[2];

	switch (code)
	{
	case E_X:
		res = vecEps[0];
		break;
	case E_Y:
		res = vecEps[1];
		break;
	case E_Z:
		res = vecEps[2];
		break;
	case E_XY:
		res = vecEps[3];
		break;
	case E_YZ:
		res = vecEps[4];
		break;
	case E_XZ:
		res = vecEps[5];
		break;
	case E_VOL:
		res = J;
		break;
	case S_X:
		res = vecSig[0];
		break;
	case S_Y:
		res = vecSig[1];
		break;
	case S_Z:
		res = vecSig[2];
		break;
	case S_XY:
		res = vecSig[3];
		break;
	case S_YZ:
		res = vecSig[4];
		break;
	case S_XZ:
		res = vecSig[5];
		break;
	case S_P:
		res =  storage->get_qi_n(-(int32)el, 0);;
		break;
	default:
		//warning("Element::getComponent: error in code %d", code);
		res=0.0;
	}
	return res;
}
//

//------------------MIXED_8N_3D_P0_NL--------------------
void MIXED_8N_3D_P0_NL::pre (uint32 el, FE_Storage_Interface *storage)
{
	S.assign(npow(n_int(),n_dim()), Vec<6>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
	C.assign(npow(n_int(),n_dim()), Vec<6>(1.0, 1.0, 1.0, 0.0, 0.0, 0.0));
	O.assign(npow(n_int(),n_dim()), Vec<9>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

	if (det.size()==0) 
	{
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->element_nodes(el, nodes_p);	
		make_Jacob(el, nodes_p);
	}
}
void MIXED_8N_3D_P0_NL::build (uint32 el, FE_Storage_Interface *storage)
{
//построение матрицы жесткости и вектора нагрузки для элемента
	Mat<24,24> Kuu; //матрица жесткости перед вектором перемещений
	double Kpp = 0.0;
	double Fp = 0.0;
	Mat<24,1> Kup;
	Vec<24> Fu; //вектор узловых сил элемента
	Vec<24> F_ext; //вектор внешних сил (пока не подсчитывается)
	
	double p_e = storage->get_qi_n(-(int32)el, 0);
	double dWt; //множитель при суммировании квадратур Гаусса
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		dWt = g_weight(nPoint);
		Mat<6,6> matE_c = storage->getMaterial().mat_E_c(ANALYSIS_3D, C[nPoint].ptr(), p_e).toMat<6,6>();
		Mat<6,1> matE_p;
		matE_p = storage->getMaterial().mat_E_p(ANALYSIS_3D, C[nPoint].ptr()).toVec<6>();
		Mat<6,24> matB_L = make_B_L(nPoint);
		Mat<9,9> matS = make_S(nPoint); //матрица S для матричного умножения
		Mat<6,9> matO = make_Omega(nPoint);
		Mat<9,24> matB_NL = make_B_NL(nPoint);

		Mat<6,24> matB_L1 = matO * matB_NL * 2;
		Mat<6,24> matB = matB_L + matB_L1; //складываем матрицы распределения лин. деф.
		Kuu += (matB.transpose() * matE_c * matB * 0.5 + matB_NL.transpose() * matS * matB_NL)*dWt;
		//double J = (1.0+O[nPoint][0])*((1.0+O[nPoint][4])*(1.0+O[nPoint][8])-O[nPoint][5]*O[nPoint][7]) - O[nPoint][1]*(O[nPoint][3]*(1.0+O[nPoint][8])-O[nPoint][6]*O[nPoint][5])+O[nPoint][2]*(O[nPoint][3]*O[nPoint][7]-O[nPoint][6]*(1.0+O[nPoint][4]));
		double J = storage->getMaterial().getJ(ANALYSIS_3D, C[nPoint].ptr());
		Fp += -(J - 1 - p_e/storage->getMaterial().getK0())*dWt;
		Fu += (matB.transpose() * S[nPoint] * (dWt*(-0.5)));
		Kup+= matB.transpose()*matE_p*(dWt*0.5);
		Kpp += -1.0/storage->getMaterial().getK0()*dWt;
	}//прошлись по всем точкам интегрирования
	

	Mat<25,25> Ke;
	Vec<25> F_e; //вектор правых частей элемента
	//сборка в одну матрицу
	for (uint16 i=0; i < 24; i++)
		for (uint16 j=0; j < 24; j++)
			Ke[i][j] = Kuu[i][j];
	for (uint16 i=0; i<24; i++)
		Ke[i][24] = Kup[i][0];
	for (uint16 i=0; i<24; i++)
		Ke[24][i] = Kup[i][0];
	Ke[24][24] = Kpp;
	for (uint16 i=0; i < 24; i++)
		F_e[i] = Fu[i];
	F_e[24] = Fp;

	assemble(el,Ke, F_e, storage);
}
//
void MIXED_8N_3D_P0_NL::update (uint32 el, FE_Storage_Interface *storage)
{
	Vec<25> Un; //вектор решений для степеней свобод элемента и его узлов
	// получаем вектор перемещений элемента из общего решения
	storage->get_q_e(el, Un.ptr());
	Vec<24> U;
	for (uint16 i=0;i<24;i++)
		U[i] = Un[i];
	double p_e = Un[24];

	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		O[nPoint] = make_B_NL(nPoint) * U;
		C[nPoint][0] = 1.0+2*O[nPoint][0]+pow(O[nPoint][0],2)+pow(O[nPoint][3],2)+pow(O[nPoint][6],2);	//C11
		C[nPoint][1] = 1.0+2*O[nPoint][4]+pow(O[nPoint][1],2)+pow(O[nPoint][4],2)+pow(O[nPoint][7],2);	//C22
		C[nPoint][2] = 1.0+2*O[nPoint][8]+pow(O[nPoint][2],2)+pow(O[nPoint][5],2)+pow(O[nPoint][8],2);	//C33
		C[nPoint][3] = O[nPoint][1]+O[nPoint][3]+O[nPoint][0]*O[nPoint][1]+O[nPoint][3]*O[nPoint][4]+O[nPoint][6]*O[nPoint][7];  //C12
		C[nPoint][4] = O[nPoint][5]+O[nPoint][7]+O[nPoint][1]*O[nPoint][2]+O[nPoint][4]*O[nPoint][5]+O[nPoint][7]*O[nPoint][8];	//C23
		C[nPoint][5] = O[nPoint][2]+O[nPoint][6]+O[nPoint][0]*O[nPoint][2]+O[nPoint][3]*O[nPoint][5]+O[nPoint][6]*O[nPoint][8];	//C13
		//восстановление напряжений Пиолы-Кирхгоффа из текущего состояния
		S[nPoint] = storage->getMaterial().getS(ANALYSIS_3D, C[nPoint].ptr(), p_e).toVec<6>();
	}
}

Mat<6,24> MIXED_8N_3D_P0_NL::make_B_L (uint16 nPoint)
{
	Mat<6,24> B_L;
	for (uint16 i=0; i < 8; i++)
	{
		B_L[0][i*3+0] = 2*NjXi[nPoint][0][i];
		B_L[1][i*3+1] = 2*NjXi[nPoint][1][i];
		B_L[2][i*3+2] = 2*NjXi[nPoint][2][i];
		B_L[3][i*3+0] = 2*NjXi[nPoint][1][i];
		B_L[3][i*3+1] = 2*NjXi[nPoint][0][i];
		B_L[4][i*3+1] = 2*NjXi[nPoint][2][i];
		B_L[4][i*3+2] = 2*NjXi[nPoint][1][i];
		B_L[5][i*3+0] = 2*NjXi[nPoint][2][i];
		B_L[5][i*3+2] = 2*NjXi[nPoint][0][i];
	}
	return B_L;
}

Mat<9,24>  MIXED_8N_3D_P0_NL::make_B_NL (uint16 nPoint)
{
	Mat<9,24> B_NL;
	for (uint16 i=0; i < 8; i++)
	{
		B_NL[0][i*3+0] = NjXi[nPoint][0][i];
		B_NL[1][i*3+0] = NjXi[nPoint][1][i];
		B_NL[2][i*3+0] = NjXi[nPoint][2][i];
		B_NL[3][i*3+1] = NjXi[nPoint][0][i];
		B_NL[4][i*3+1] = NjXi[nPoint][1][i];
		B_NL[5][i*3+1] = NjXi[nPoint][2][i];
		B_NL[6][i*3+2] = NjXi[nPoint][0][i];
		B_NL[7][i*3+2] = NjXi[nPoint][1][i];
		B_NL[8][i*3+2] = NjXi[nPoint][2][i];
	}
	return B_NL;
}

Mat<9,9> MIXED_8N_3D_P0_NL::make_S (uint16 nPoint)
{
	Mat<9,9> matS;
	for (uint16 i=0; i < 3; i++)
	{
		matS[i*3+0][i*3+0] = S[nPoint][0];
		matS[i*3+0][i*3+1] = S[nPoint][3];
		matS[i*3+0][i*3+2] = S[nPoint][5];
		matS[i*3+1][i*3+0] = S[nPoint][3];
		matS[i*3+1][i*3+1] = S[nPoint][1];
		matS[i*3+1][i*3+2] = S[nPoint][4];
		matS[i*3+2][i*3+0] = S[nPoint][5];
		matS[i*3+2][i*3+1] = S[nPoint][4];
		matS[i*3+2][i*3+2] = S[nPoint][2];
	}
	return matS;
}

Mat<6,9> MIXED_8N_3D_P0_NL::make_Omega (uint16 nPoint)
{
	Mat<6,9> Omega;
	for (uint16 i=0; i < 3; i++)
	{
		Omega[0][i*3+0] = O[nPoint][0+i*3];
		Omega[1][i*3+1] = O[nPoint][1+i*3];
		Omega[2][i*3+2] = O[nPoint][2+i*3];
		Omega[3][i*3+0] = O[nPoint][1+i*3];
		Omega[3][i*3+1] = O[nPoint][0+i*3];
		Omega[4][i*3+1] = O[nPoint][2+i*3];
		Omega[4][i*3+2] = O[nPoint][1+i*3];
		Omega[5][i*3+0] = O[nPoint][2+i*3];
		Omega[5][i*3+2] = O[nPoint][0+i*3];
	}
	return Omega;
}

double MIXED_8N_3D_P0_NL::getComponent(uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage)
{
	//////////////////////////////////////
	//see codes in sys.h
	//gp - needed gauss point 

	assert (gp < n_int()*n_int()*n_int());
	double res = 0.0;

	if (code >= E_X && code <= E_3)
	{
		Mat<3,3> matC(C[gp][0],C[gp][3],C[gp][5],
				  C[gp][3],C[gp][1],C[gp][4],
				  C[gp][5],C[gp][4],C[gp][2]);
		 
		switch (code)
		{
		case E_X:
			res = (matC[0][0]-1)*0.5;
			break;
		case E_Y:
			res = (matC[1][1]-1)*0.5;
			break;
		case E_Z:
			res = (matC[2][2]-1)*0.5;
			break;
		case E_XY:
			res = matC[0][1]*0.5;
			break;
		case E_YZ:
			res = matC[1][2]*0.5;
			break;
		case E_XZ:
			res = matC[0][2]*0.5;
			break;
		case E_VOL:
			res = sqrt(matC.det());
			break;
		case E_1:
			res = sqrt(matC.eigenvalues()[0]);
			break;
		case E_2:
			res = sqrt(matC.eigenvalues()[1]);
			break;
		case E_3:
			res = sqrt(matC.eigenvalues()[2]);
			break;
		}
	}

	if (code >= S_X && code <= S_3)
	{
		Mat<3,3> matX;
		Mat<3,3> matS;
		Mat<3,3> matT;

		matX[0][0] = 1+O[gp][0];
		matX[0][1] = O[gp][1];
		matX[0][2] = O[gp][2];
		matX[1][0] = O[gp][3];
		matX[1][1] = 1+O[gp][4];
		matX[1][2] = O[gp][5];
		matX[2][0] = O[gp][6];
		matX[2][1] = O[gp][7];
		matX[2][2] = 1+O[gp][8];

		matS[0][0] = S[gp][0];
		matS[0][1] = S[gp][3];
		matS[0][2] = S[gp][5];
		matS[1][0] = S[gp][3];
		matS[1][1] = S[gp][1];
		matS[1][2] = S[gp][4];
		matS[2][0] = S[gp][5];
		matS[2][1] = S[gp][4];
		matS[2][2] = S[gp][2];

		matT = matX * matS * matX.transpose() * (1.0f/matX.det());

		switch (code)
		{
		case S_X:
			res = matT[0][0];
			break;
		case S_Y:
			res = matT[1][1];
			break;
		case S_Z:
			res = matT[2][2];
			break;
		case S_XY:
			res = matT[0][1];
			break;
		case S_YZ:
			res = matT[1][2];
			break;
		case S_XZ:
			res = matT[0][2];
			break;
		case S_P:
			res = storage->get_qi_n(-(int32)el, 0);
			break;
		case S_1:
			res = matT.eigenvalues()[0];
			break;
		case S_2:
			res = matT.eigenvalues()[1];
			break;
		case S_3:
			res = matT.eigenvalues()[2];
			break;
		}
	}
	return res;
}




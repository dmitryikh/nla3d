#include "element_MIXED_8N_3D_P0.h"

const mat_comp MIXED_8N_3D_P0::components[6] = {M_XX, M_YY, M_ZZ, M_XY, M_YZ, M_XZ};
const uint16 MIXED_8N_3D_P0::num_components = 6;


void MIXED_8N_3D_P0::pre (uint32 el, FE_Storage_Interface *storage)
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


void MIXED_8N_3D_P0::build (uint32 el, FE_Storage_Interface *storage)
{
//построение матрицы жесткости и вектора нагрузки для элемента
	double Kpp = 0.0;
	double Fp = 0.0;
	
	Vec<24> Kup;
	Vec<24> Fu; //вектор узловых сил элемента
	Vec<24> F_ext; //вектор внешних сил (пока не подсчитывается)
	Material *mat = (Material*) GlobStates.getptr(States::PTR_CURMATER);
	double k = mat->getK0();
	MatSym<6> matD_d;
	Vec<6> vecD_p;
	Vec<6> vecC;

	Mat2<6,24> matB;
	MatSym<9> matS;
	Mat2<6,9> matO;
	Mat2<9,24> matB_NL;
	MatSym<24> Kuu; //матрица жесткости перед вектором перемещений
	double p_e = storage->get_qi_n(-(int32)el, 0);
	GlobStates.setdouble(States::DOUBLE_HYDPRES, p_e);
	double dWt; //множитель при суммировании квадратур Гаусса

	Kuu.zeros();
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		dWt = g_weight(nPoint);

		//C[0] - C11	C[1] - C22	C[2] - C33	C[3] - C12	C[4] - C23	C[5] - C13
		vecC[M_XX] = C[nPoint][0];
		vecC[M_YY] = C[nPoint][1];
		vecC[M_ZZ] = C[nPoint][2];
		vecC[M_XY] = C[nPoint][3];
		vecC[M_YZ] = C[nPoint][4];
		vecC[M_XZ] = C[nPoint][5];
		mat->getDdDp_UP(num_components, components, vecC.ptr(), matD_d.ptr(), vecD_p.ptr());
		double J = Material::getJ(vecC.ptr());
		//Mat<6,6> matE_c = storage->getMaterial().mat_E_c(ANALYSIS_3D, C[nPoint].ptr(), p_e).toMat<6,6>();
		//Mat<6,1> matE_p;
		//matE_p = storage->getMaterial().mat_E_p(ANALYSIS_3D, C[nPoint].ptr()).toVec<6>();
		matB.zeros();
		matS.zeros();
		matO.zeros();
		matB_NL.zeros();

		make_B_L(nPoint, matB);
		make_S(nPoint, matS); //матрица S для матричного умножения
		make_Omega(nPoint, matO);
		make_B_NL(nPoint, matB_NL);
		matABprod(matO, matB_NL, 2.0, matB);
		//Mat<6,9> matO = make_Omega(nPoint);
		//Mat<9,24> matB_NL = make_B_NL(nPoint);
		//was: Mat<6,24> matB_L1 = matO * matB_NL * 2;
		//Mat2<6,24> matB_L1;
		//was: Mat<6,24> matB = matB_L + matB_L1; //складываем матрицы распределения лин. деф.
		//was: Kuu += (matB.transpose() * matE_c * matB * 0.5 + matB_NL.transpose() * matS * matB_NL)*dWt;
		matBTDBprod(matB, matD_d, 0.5*dWt, Kuu);
		matBTDBprod(matB_NL, matS, dWt, Kuu);

		//was: Fu += (matB.transpose() * S[nPoint] * (dWt*(-0.5)));
		matBTVprod(matB,S[nPoint], -0.5*dWt, Fu);

		//was: Kup+= matB.transpose()*matE_p*(dWt*0.5);
		matBTVprod(matB,vecD_p, 0.5*dWt, Kup);

		Fp += -(J - 1 - p_e/k)*dWt;
		Kpp += -1.0/k*dWt;
	}//прошлись по всем точкам интегрирования
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);
	assemble3(el,Kuu, Kup, Kpp, Fu,Fp,  storage);

	//Mat<25,25> Ke;
	//Vec<25> F_e; //вектор правых частей элемента
	////сборка в одну матрицу
	//for (uint16 i=0; i < 24; i++)
	//	for (uint16 j=0; j < 24; j++)
	//		Ke[i][j] = Kuu[i][j];
	//for (uint16 i=0; i<24; i++)
	//	Ke[i][24] = Kup[i][0];
	//for (uint16 i=0; i<24; i++)
	//	Ke[24][i] = Kup[i][0];
	//Ke[24][24] = Kpp;
	//for (uint16 i=0; i < 24; i++)
	//	F_e[i] = Fu[i];
	//F_e[24] = Fp;

	//assemble(el,Ke, F_e, storage);
}


void MIXED_8N_3D_P0::update (uint32 el, FE_Storage_Interface *storage)
{
	Vec<25> Un; //вектор решений для степеней свобод элемента и его узлов
	// получаем вектор перемещений элемента из общего решения
	storage->get_q_e(el, Un.ptr());
	Vec<24> U;
	Mat2<9,24> B_NL;
	Vec<6> vecC;
	for (uint16 i=0;i<24;i++)
		U[i] = Un[i];
	double p_e = Un[24];
	GlobStates.setdouble(States::DOUBLE_HYDPRES, p_e);
	Material *mat = (Material*) GlobStates.getptr(States::PTR_CURMATER);
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		B_NL.zeros();
		make_B_NL(nPoint, B_NL);
		O[nPoint].zeros();
		matBVprod(B_NL, U, 1.0, O[nPoint]);
		//O[nPoint] = make_B_NL(nPoint) * U;
		C[nPoint][0] = 1.0+2*O[nPoint][0]+pow(O[nPoint][0],2)+pow(O[nPoint][3],2)+pow(O[nPoint][6],2);	//C11
		C[nPoint][1] = 1.0+2*O[nPoint][4]+pow(O[nPoint][1],2)+pow(O[nPoint][4],2)+pow(O[nPoint][7],2);	//C22
		C[nPoint][2] = 1.0+2*O[nPoint][8]+pow(O[nPoint][2],2)+pow(O[nPoint][5],2)+pow(O[nPoint][8],2);	//C33
		C[nPoint][3] = O[nPoint][1]+O[nPoint][3]+O[nPoint][0]*O[nPoint][1]+O[nPoint][3]*O[nPoint][4]+O[nPoint][6]*O[nPoint][7];  //C12
		C[nPoint][4] = O[nPoint][5]+O[nPoint][7]+O[nPoint][1]*O[nPoint][2]+O[nPoint][4]*O[nPoint][5]+O[nPoint][7]*O[nPoint][8];	//C23
		C[nPoint][5] = O[nPoint][2]+O[nPoint][6]+O[nPoint][0]*O[nPoint][2]+O[nPoint][3]*O[nPoint][5]+O[nPoint][6]*O[nPoint][8];	//C13
		//восстановление напряжений Пиолы-Кирхгоффа из текущего состояния

		vecC[M_XX] = C[nPoint][0];
		vecC[M_YY] = C[nPoint][1];
		vecC[M_ZZ] = C[nPoint][2];
		vecC[M_XY] = C[nPoint][3];
		vecC[M_YZ] = C[nPoint][4];
		vecC[M_XZ] = C[nPoint][5];
		mat->getS_UP (num_components, components, vecC.ptr(), S[nPoint].ptr());
	}
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);
}


void MIXED_8N_3D_P0::make_B_L (uint16 nPoint, Mat2<6,24> &B)
{
	double *B_L = B.ptr();
	for (uint16 i=0; i < 8; i++)
	{
		B_L[0*24+(i*3+0)] += 2*NjXi[nPoint][0][i];
		B_L[1*24+(i*3+1)] += 2*NjXi[nPoint][1][i];
		B_L[2*24+(i*3+2)] += 2*NjXi[nPoint][2][i];
		B_L[3*24+(i*3+0)] += 2*NjXi[nPoint][1][i];
		B_L[3*24+(i*3+1)] += 2*NjXi[nPoint][0][i];
		B_L[4*24+(i*3+1)] += 2*NjXi[nPoint][2][i];
		B_L[4*24+(i*3+2)] += 2*NjXi[nPoint][1][i];
		B_L[5*24+(i*3+0)] += 2*NjXi[nPoint][2][i];
		B_L[5*24+(i*3+2)] += 2*NjXi[nPoint][0][i];
	}
}


void MIXED_8N_3D_P0::make_B_NL (uint16 nPoint,  Mat2<9,24> &B)
{
	double *B_NL = B.ptr();
	for (uint16 i=0; i < 8; i++)
	{
		B_NL[0*24+(i*3+0)] += NjXi[nPoint][0][i];
		B_NL[1*24+(i*3+0)] += NjXi[nPoint][1][i];
		B_NL[2*24+(i*3+0)] += NjXi[nPoint][2][i];
		B_NL[3*24+(i*3+1)] += NjXi[nPoint][0][i];
		B_NL[4*24+(i*3+1)] += NjXi[nPoint][1][i];
		B_NL[5*24+(i*3+1)] += NjXi[nPoint][2][i];
		B_NL[6*24+(i*3+2)] += NjXi[nPoint][0][i];
		B_NL[7*24+(i*3+2)] += NjXi[nPoint][1][i];
		B_NL[8*24+(i*3+2)] += NjXi[nPoint][2][i];
	}
}



void MIXED_8N_3D_P0::make_S (uint16 nPoint, MatSym<9> &B)
{
	double *Sp = B.ptr();
	Sp[0]	+= S[nPoint][0];
	Sp[1]	+= S[nPoint][3];
	Sp[2]	+= S[nPoint][5];
	Sp[9]	+= S[nPoint][1];
	Sp[10]	+= S[nPoint][4];
	Sp[17]	+= S[nPoint][2];

	Sp[24]	+= S[nPoint][0];
	Sp[25]	+= S[nPoint][3];
	Sp[26]	+= S[nPoint][5];
	Sp[30]	+= S[nPoint][1];
	Sp[31]	+= S[nPoint][4];
	Sp[35]	+= S[nPoint][2];

	Sp[39]	+= S[nPoint][0];
	Sp[40]	+= S[nPoint][3];
	Sp[41]	+= S[nPoint][5];
	Sp[42]	+= S[nPoint][1];
	Sp[43]	+= S[nPoint][4];
	Sp[44]	+= S[nPoint][2];

	/*Mat<9,9> matS;
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
	return matS;*/
}


void MIXED_8N_3D_P0::make_Omega (uint16 nPoint, Mat2<6,9> &B)
{
	double *Omega = B.ptr();

	for (uint16 i=0; i < 3; i++)
	{
		Omega[0*9+(i*3+0)] = O[nPoint][0+i*3];
		Omega[1*9+(i*3+1)] = O[nPoint][1+i*3];
		Omega[2*9+(i*3+2)] = O[nPoint][2+i*3];
		Omega[3*9+(i*3+0)] = O[nPoint][1+i*3];
		Omega[3*9+(i*3+1)] = O[nPoint][0+i*3];
		Omega[4*9+(i*3+1)] = O[nPoint][2+i*3];
		Omega[4*9+(i*3+2)] = O[nPoint][1+i*3];
		Omega[5*9+(i*3+0)] = O[nPoint][2+i*3];
		Omega[5*9+(i*3+2)] = O[nPoint][0+i*3];
	}
}


double MIXED_8N_3D_P0::getComponent(uint16 gp, el_component code, uint32 el, FE_Storage_Interface *storage)
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


Mat<3,3> MIXED_8N_3D_P0::getTensor (uint16 gp, el_tensor code, uint32 el, FE_Storage_Interface *storage) {
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

	Mat<3,3> matC;
	matC = Q.transpose()*matT*Q;



	switch (code)
	{
	case TENS_COUCHY:
		return matC;
		break;
	}

	return Mat<3,3>(0.0,0.0,0.0,
					0.0,0.0,0.0,
					0.0,0.0,0.0);
}

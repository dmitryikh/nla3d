#include "elements/SOLID80.h"

const tensorComponents ElementSOLID80::components[6] = {M_XX, M_YY, M_ZZ, M_XY, M_YZ, M_XZ};
const uint16 ElementSOLID80::num_components = 6;


void ElementSOLID80::pre()
{
	S.assign(npow(n_int(),n_dim()), Vec<6>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
	C.assign(npow(n_int(),n_dim()), Vec<6>(1.0, 1.0, 1.0, 0.0, 0.0, 0.0));
	O.assign(npow(n_int(),n_dim()), Vec<9>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

	if (det.size()==0) {
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->element_nodes(getElNum(), nodes_p);	
		make_Jacob(getElNum(), nodes_p);
	}
}


void ElementSOLID80::build()
{
//построение матрицы жесткости и вектора нагрузки для элемента
	double Kpp = 0.0;
	double Fp = 0.0;
	
	Vec<24> Kup;
	Vec<24> Fu; //вектор узловых сил элемента
	Vec<24> F_ext; //вектор внешних сил (пока не подсчитывается)
  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> ((Material*)GlobStates.getptr(States::PTR_CURMATER));
  if (mat == NULL) {
    error("SOLLID81::build: material is not derived from Mat_Hyper_Isotrop_General");
  }
	double k = mat->getK();
	MatSym<6> matD_d;
	Vec<6> vecD_p;
	Vec<6> vecC;

	Mat2<6,24> matB;
	MatSym<9> matS;
	Mat2<6,9> matO;
	Mat2<9,24> matB_NL;
	MatSym<24> Kuu; //матрица жесткости перед вектором перемещений
	double p_e = storage->get_qi_n(-(int32)getElNum(), 0);
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
	assemble3(Kuu, Kup, Kpp, Fu,Fp);
}


void ElementSOLID80::update()
{
	Vec<25> Un; //вектор решений для степеней свобод элемента и его узлов
	// получаем вектор перемещений элемента из общего решения
	storage->get_q_e(getElNum(), Un.ptr());
	Vec<24> U;
	Mat2<9,24> B_NL;
	Vec<6> vecC;
	for (uint16 i=0;i<24;i++)
		U[i] = Un[i];
	double p_e = Un[24];
	GlobStates.setdouble(States::DOUBLE_HYDPRES, p_e);
  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> ((Material*)GlobStates.getptr(States::PTR_CURMATER));
  if (mat == NULL) {
    error("ElementSOLID80::update: material is not derived from Mat_Hyper_Isotrop_General");
  }
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


void ElementSOLID80::make_B_L (uint16 nPoint, Mat2<6,24> &B)
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


void ElementSOLID80::make_B_NL (uint16 nPoint,  Mat2<9,24> &B)
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



void ElementSOLID80::make_S (uint16 nPoint, MatSym<9> &B)
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


void ElementSOLID80::make_Omega (uint16 nPoint, Mat2<6,9> &B)
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

void ElementSOLID80::getScalar(double& scalar, el_component code, uint16 gp, const double scale) {
	//see codes in sys.h
	//gp - needed gauss point 
  if (gp == GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 nPoint = 0; nPoint < npow(n_int(),n_dim()); nPoint ++) {
      dWt = g_weight(nPoint);
      getScalar(scalar, code, nPoint, dWt/dWtSum*scale );
    }
    return;
  }
	assert (npow(n_int(),n_dim()));
  switch (code) {
		case S_P:
			scalar += storage->get_qi_n(-(int32)getElNum(), 0) * scale;
			break;
    default:
      error("ElementSOLID80::getScalar: no data for code %d", code);
  }
}

//return a tensor in a global coordinate system
void  ElementSOLID80::getTensor(MatSym<3>& tensor, el_tensor code, uint16 gp, const double scale) {
  if (gp == GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 nPoint = 0; nPoint < npow(n_int(),n_dim()); nPoint ++) {
      dWt = g_weight(nPoint);
      getTensor(tensor, code, nPoint, dWt/dWtSum*scale);
    }
    return;
  }
	assert (npow(n_int(),n_dim()));

  Mat2<3,3> matF;
  MatSym<3> matS;
  double J;

	switch (code) {
    case TENS_COUCHY:

      //matF^T  
      matF.data[0][0] = 1+O[gp][0];
      matF.data[1][0] = O[gp][1];
      matF.data[2][0] = O[gp][2];
      matF.data[0][1] = O[gp][3];
      matF.data[1][1] = 1+O[gp][4];
      matF.data[2][1] = O[gp][5];
      matF.data[0][2] = O[gp][6];
      matF.data[1][2] = O[gp][7];
      matF.data[2][2] = 1+O[gp][8];

      J = matF.data[0][0]*(matF.data[1][1]*matF.data[2][2]-matF.data[1][2]*matF.data[2][1])-matF.data[0][1]*(matF.data[1][0]*matF.data[2][2]-matF.data[1][2]*matF.data[2][0])+matF.data[0][2]*(matF.data[1][0]*matF.data[2][1]-matF.data[1][1]*matF.data[2][0]);

      matS.data[0] = S[gp][0];
      matS.data[1] = S[gp][3];
      matS.data[2] = S[gp][5];
      matS.data[3] = S[gp][1];
      matS.data[4] = S[gp][4];
      matS.data[5] = S[gp][2];
      matBTDBprod (matF, matS, 1.0/J*scale, tensor); //Symmetric Couchy tensor
      break;
    case TENS_PK2:
      tensor.data[0] += S[gp][0]*scale;
      tensor.data[1] += S[gp][3]*scale;
      tensor.data[2] += S[gp][5]*scale;
      tensor.data[3] += S[gp][1]*scale;
      tensor.data[4] += S[gp][4]*scale;
      tensor.data[5] += S[gp][2]*scale;
      break;
    case TENS_C:
      tensor.data[0] += C[gp][0]*scale;
      tensor.data[1] += C[gp][3]*scale;
      tensor.data[2] += C[gp][5]*scale;
      tensor.data[3] += C[gp][1]*scale;
      tensor.data[4] += C[gp][4]*scale;
      tensor.data[5] += C[gp][2]*scale;
      break;
    case TENS_E:
      tensor.data[0] += (C[gp][0]-1.0)*0.5*scale;
      tensor.data[1] += C[gp][3]*0.5*scale;
      tensor.data[2] += C[gp][5]*0.5*scale;
      tensor.data[3] += (C[gp][1]-1.0)*0.5*scale;
      tensor.data[4] += C[gp][4]*0.5*scale;
      tensor.data[5] += (C[gp][2]-1.0)*0.5*scale;
      break;
    default:
      error("ElementSOLID80::getTensor: no data for code %d", code);
	}
}

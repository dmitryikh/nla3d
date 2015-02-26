#include "element_MIXED_4N_2D_P0.h"

const mat_comp MIXED_4N_2D_P0::components[3] = {M_XX, M_YY, M_XY};
const uint16 MIXED_4N_2D_P0::num_components = 3;

//------------------MIXED_4N_2D_P0--------------------
void MIXED_4N_2D_P0::pre ()
{
	S.assign(npow(n_int(),n_dim()), Vec<3>(0.0f, 0.0f, 0.0f));
	C.assign(npow(n_int(),n_dim()), Vec<3>(1.0f, 1.0f, 0.0f));
	O.assign(npow(n_int(),n_dim()), Vec<4>(0.0f, 0.0f, 0.0f, 0.0f));
	if (det.size()==0) 
	{
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->element_nodes(getElNum(), nodes_p);	
		make_Jacob(getElNum(), nodes_p);
	}
}
void MIXED_4N_2D_P0::build ()
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
	double p_e = storage->get_qi_n(-(int32)getElNum(), 0);
	GlobStates.setdouble(States::DOUBLE_HYDPRES, p_e);
  Material *mat = (Material*) GlobStates.getptr(States::PTR_CURMATER);
	double k = mat->getK0();
	double dWt; //Gaussian quadrature
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		dWt = g_weight(nPoint);
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
		Fp += (J - 1 - p_e/k)*dWt;
		Qe += (matB.transpose() * S[nPoint] * dWt);
		Kup+= matB.transpose()*matE_p *dWt;
		Kpp -= 1.0/k*dWt;

	}// loop over intergration points
	
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
	assemble(getElNum(), Ke, Fe, storage);
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
void MIXED_4N_2D_P0::update ()
{
	Vec<9> U; //вектор перемещений элемента
	// получаем вектор перемещений элемента из общего решения
	storage->get_q_e(getElNum(), U.ptr());
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
    //all meterial functions are waiting [C] for 3D case. So we need to use CVec here.
    CVec[M_XX] = C[nPoint][0];
    CVec[M_YY] = C[nPoint][1];
    CVec[M_XY] = C[nPoint][2];
		mat->getS_UP (num_components, components, CVec.ptr(), S[nPoint].ptr());
	}
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);
}

void MIXED_4N_2D_P0::getScalar(double& scalar, el_component code, uint16 gp, const double scale) {
	//see codes in sys.h
	//gp - needed gauss point 
  if {gp == GP_MEAN} { //need to average result over the element
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
      error("MIXED_4N_2D_P0::getScalar: no data for code %d", code);
  }
}

//return a tensor in a global coordinate system
void  MIXED_4N_2D_P0::getTensor(MatSym<3>& tensor, el_tensor code, uint16 gp, const double scale) {
  if {gp == GP_MEAN} { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 nPoint = 0; nPoint < npow(n_int(),n_dim()); nPoint ++) {
      dWt = g_weight(nPoint);
      getTensor(tensor, code, nPoint, dWt/dWtSum*scale);
    }
    return;
  }
	assert (npow(n_int(),n_dim()));

  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
  CVec[M_XX] = C[gp][0];
  CVec[M_YY] = C[gp][1];
  CVec[M_XY] = C[gp][2];
  //was
  //matS[0][0] = S[gp][0]; //SX
  //matS[0][1] = S[gp][2]; //SXY
  //matS[0][2] = 0.0;
  //matS[1][0] = S[gp][2]; //SXY
  //matS[1][1] = S[gp][1]; //SYY
  //matS[1][2] = 0.0;
  //matS[2][0] = 0.0;
  //matS[2][1] = 0.0;
  //matS[2][2] = 0.0;//it's worng. Actually we have SZZ
	//matT = matX * matS * matX.transpose() * (1.0f/matX.det());

	switch (code) {
    case TENS_COUCHY:
      Mat2<3,3> matF;
      MatSym<3> matS;
      double J;

      matF.zero();
      matF[0][0] = 1+O[gp][0];
      matF[0][1] = O[gp][1];
      matF[1][0] = O[gp][2];
      matF[1][1] = 1+O[gp][3];
      matF[2][2] = 1;

      J = matF.data[0][0]*(matF.data[1][1]*matF.data[2][2]-matF.data[1][2]*matF.data[2][1])-matF.data[0][1]*(matF.data[1][0]*matF.data[2][2]-matF.data[1][2]*matF.data[2][0])+matF.data[0][2]*(matF.data[1][0]*matF.data[2][1]-matF.data[1][1]*matF.data[2][0]);
      //In order to complete matS (3x3 symmetric matrix, PK2 tensor) we need 
      //to know S33 component: 
      //1) One solution is to calculate S33 on every solution step
      //and store it in S[nPoint] vector.
      //2) Second solution is to resotre S33 right here.
      //Now 2) is working.
      storage->getMaterial().getS_UP (6, MatCompsGlobal, CVec.ptr(), matS.data);
      matBTDBprod (matF, matS, 1.0/J, tensor); //Symmetric Couchy tensor
      break;
    case TENS_PK2:
      MatSym<3> matS;
      storage->getMaterial().getS_UP (6, MatCompsGlobal, CVec.ptr(), matS.data);
      tensor.data[0] += matS.data[0]*scale;
      tensor.data[1] += matS.data[1]*scale;
      tensor.data[2] += matS.data[2]*scale;
      tensor.data[3] += matS.data[3]*scale;
      tensor.data[4] += matS.data[4]*scale;
      tensor.data[5] += matS.data[5]*scale;
      break;
    case TENS_C:
      tensor.data[0] += CVec[M_XX]*scale;
      tensor.data[1] += CVec[M_XY]*scale;
      tensor.data[2] += CVec[M_XZ]*scale;
      tensor.data[3] += CVec[M_YY]*scale;
      tensor.data[4] += CVec[M_YZ]*scale;
      tensor.data[5] += CVec[M_ZZ]*scale;
      break;
    case TENS_E:
      tensor.data[0] += (CVec[M_XX]-1.0)*0.5*scale;
      tensor.data[1] += CVec[M_XY]*0.5*scale;
      tensor.data[2] += CVec[M_XZ]*0.5*scale;
      tensor.data[3] += (CVec[M_YY]-1.0)*0.5*scale;
      tensor.data[4] += CVec[M_YZ]*0.5*scale;
      tensor.data[5] += (CVec[M_ZZ]-1.0)*0.5*scale;
      break;
    default:
      error("MIXED_4N_2D_P0::getTensor: no data for code %d", code);
	}
}

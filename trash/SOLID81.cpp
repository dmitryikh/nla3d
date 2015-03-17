#include "elements/SOLID81.h"

const tensorComponents ElementSOLID81::components[6] = {M_XX, M_XY, M_XZ, M_YY, M_YZ, M_ZZ};
const uint16 ElementSOLID81::num_components = 6;


void ElementSOLID81::pre()
{
	S.assign(npow(n_int(),n_dim()), Vec<6>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
	C.assign(npow(n_int(),n_dim()), Vec<6>(1.0, 0.0, 0.0, 1.0, 0.0, 1.0));
	O.assign(npow(n_int(),n_dim()), Vec<9>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

	if (det.size()==0) {
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->element_nodes(getElNum(), nodes_p);	
		make_Jacob(getElNum(), nodes_p);
	}
}


void ElementSOLID81::build()
{
//построение матрицы жесткости и вектора нагрузки для элемента
	double Kpp = 0.0;
	double Fp = 0.0;
	
  //Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, 25, 25> Ke(25, 25);
  //Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::RowMajor, 25, 1> Fe(25, 1);
  Eigen::MatrixXd Ke(25, 25);
  Eigen::MatrixXd Fe(25, 1);

  Eigen::Matrix<double, 6, 6> matD_d;
  Eigen::Matrix<double, 6, 1> vecD_p;

  Eigen::Matrix<double, 6,24> matB;
  Eigen::Matrix<double, 9, 9> matS; //this is symmetric matrix
  Eigen::Matrix<double, 6, 9> matO;
  Eigen::Matrix<double, 9,24> matB_NL;

  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> ((Material*)GlobStates.getptr(States::PTR_CURMATER));
  if (mat == NULL) {
    error("SOLLID81::build: material is not derived from Mat_Hyper_Isotrop_General");
  }
	double k = mat->getK();
	double p_e = storage->get_qi_n(-(int32)getElNum(), 0);
	GlobStates.setdouble(States::DOUBLE_HYDPRES, p_e);
	double dWt; //множитель при суммировании квадратур Гаусса
  Ke.setZero();
  Fe.setZero();
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		dWt = g_weight(nPoint);

		mat->getDdDp_UP(6, defaultTensorComponents, C[nPoint].ptr(), matD_d, vecD_p);
		double J = Material::getJ(C[nPoint].ptr());
		matB.setZero();
		matS.setZero();
		matO.setZero();
		matB_NL.setZero();

		make_B_L(nPoint, matB);
		make_S(nPoint, matS); //matrix 9x9 with 3x3 stress tenros blocks
		make_Omega(nPoint, matO);
		make_B_NL(nPoint, matB_NL);
    matB += matO * matB_NL * 2.0;
    Ke.topLeftCorner<24,24>().triangularView<Eigen::Upper>() += 
    //Ke.topLeftCorner<24,24>() += 
      matB.transpose() * matD_d.selfadjointView<Eigen::Upper>() * matB * (0.5*dWt) +
      matB_NL.transpose() * matS.selfadjointView<Eigen::Upper>() * matB_NL * dWt;
    Fe.topLeftCorner<24,1>() += matB.transpose() * Eigen::Map<Eigen::VectorXd> (S[nPoint].ptr(), 6, 1) * (-0.5*dWt);
		Ke.topRightCorner<24,1>() += matB.transpose() *vecD_p * (0.5*dWt);

		Fe(24) += -(J - 1 - p_e/k)*dWt;
		Ke(24,24) += -1.0/k*dWt;
	}//прошлись по всем точкам интегрирования
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);
  //cout << "Ke = " << Ke << endl;
  //cout << "Fe = " << Fe << endl;
	assemble4(Ke, Fe);
}


void ElementSOLID81::update()
{
	//Vec<25> Un; //вектор решений для степеней свобод элемента и его узлов
	// получаем вектор перемещений элемента из общего решения
  Eigen::VectorXd U(25);
	storage->get_q_e(getElNum(), U.data());
  Eigen::Matrix<double, 9,24> B_NL;
	double p_e = U(24);
	GlobStates.setdouble(States::DOUBLE_HYDPRES, p_e);
  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> ((Material*)GlobStates.getptr(States::PTR_CURMATER));
  if (mat == NULL) {
    error("ElementSOLID81::update: material is not derived from Mat_Hyper_Isotrop_General");
  }
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		B_NL.setZero();
		make_B_NL(nPoint, B_NL);
		//matBVprod(B_NL, U, 1.0, O[nPoint]);
    Eigen::VectorXd tmp = B_NL*U.topLeftCorner<24,1>();
    for (uint16 i = 0; i < tmp.rows(); i++)
      O[nPoint][i] = tmp(i);
		C[nPoint][M_XX] = 1.0+2*O[nPoint][0]+pow(O[nPoint][0],2)+pow(O[nPoint][3],2)+pow(O[nPoint][6],2);	//C11
		C[nPoint][M_YY] = 1.0+2*O[nPoint][4]+pow(O[nPoint][1],2)+pow(O[nPoint][4],2)+pow(O[nPoint][7],2);	//C22
		C[nPoint][M_ZZ] = 1.0+2*O[nPoint][8]+pow(O[nPoint][2],2)+pow(O[nPoint][5],2)+pow(O[nPoint][8],2);	//C33
		C[nPoint][M_XY] = O[nPoint][1]+O[nPoint][3]+O[nPoint][0]*O[nPoint][1]+O[nPoint][3]*O[nPoint][4]+O[nPoint][6]*O[nPoint][7];  //C12
		C[nPoint][M_YZ] = O[nPoint][5]+O[nPoint][7]+O[nPoint][1]*O[nPoint][2]+O[nPoint][4]*O[nPoint][5]+O[nPoint][7]*O[nPoint][8];	//C23
		C[nPoint][M_XZ] = O[nPoint][2]+O[nPoint][6]+O[nPoint][0]*O[nPoint][2]+O[nPoint][3]*O[nPoint][5]+O[nPoint][6]*O[nPoint][8];	//C13
		//восстановление напряжений Пиолы-Кирхгоффа из текущего состояния
		mat->getS_UP (6, defaultTensorComponents, C[nPoint].ptr(), S[nPoint].ptr());
	}
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);
}


void ElementSOLID81::make_B_L (uint16 nPoint, Eigen::Ref<Eigen::MatrixXd> B) {
//was Mat2<6,24> &B
	for (uint16 i=0; i < 8; i++) {
		B(0,i*3+0) += 2*NjXi[nPoint][0][i];//exx
		B(1,i*3+0) += 2*NjXi[nPoint][1][i];//exy
		B(1,i*3+1) += 2*NjXi[nPoint][0][i];//exy
		B(2,i*3+0) += 2*NjXi[nPoint][2][i];//exz
		B(2,i*3+2) += 2*NjXi[nPoint][0][i];//exz
		B(3,i*3+1) += 2*NjXi[nPoint][1][i];//eyy
		B(4,i*3+1) += 2*NjXi[nPoint][2][i];//eyz
		B(4,i*3+2) += 2*NjXi[nPoint][1][i];//eyz
		B(5,i*3+2) += 2*NjXi[nPoint][2][i];//ezz
	}
}


void ElementSOLID81::make_B_NL (uint16 nPoint, Eigen::Ref<Eigen::MatrixXd> B) {
	for (uint16 i=0; i < 8; i++) {
		B(0,i*3+0) += NjXi[nPoint][0][i];
		B(1,i*3+0) += NjXi[nPoint][1][i];
		B(2,i*3+0) += NjXi[nPoint][2][i];
		B(3,i*3+1) += NjXi[nPoint][0][i];
		B(4,i*3+1) += NjXi[nPoint][1][i];
		B(5,i*3+1) += NjXi[nPoint][2][i];
		B(6,i*3+2) += NjXi[nPoint][0][i];
		B(7,i*3+2) += NjXi[nPoint][1][i];
		B(8,i*3+2) += NjXi[nPoint][2][i];
	}
}



void ElementSOLID81::make_S (uint16 nPoint, Eigen::Ref<Eigen::MatrixXd> SMat)
{
  SMat(0,0) += S[nPoint][M_XX];
  SMat(0,1) += S[nPoint][M_XY];
  SMat(0,2) += S[nPoint][M_XZ];
  SMat(1,1) += S[nPoint][M_YY];
  SMat(1,2) += S[nPoint][M_YZ];
  SMat(2,2) += S[nPoint][M_ZZ];

  SMat(3+0,3+0) += S[nPoint][M_XX];
  SMat(3+0,3+1) += S[nPoint][M_XY];
  SMat(3+0,3+2) += S[nPoint][M_XZ];
  SMat(3+1,3+1) += S[nPoint][M_YY];
  SMat(3+1,3+2) += S[nPoint][M_YZ];
  SMat(3+2,3+2) += S[nPoint][M_ZZ];

  SMat(6+0,6+0) += S[nPoint][M_XX];
  SMat(6+0,6+1) += S[nPoint][M_XY];
  SMat(6+0,6+2) += S[nPoint][M_XZ];
  SMat(6+1,6+1) += S[nPoint][M_YY];
  SMat(6+1,6+2) += S[nPoint][M_YZ];
  SMat(6+2,6+2) += S[nPoint][M_ZZ];
}


void ElementSOLID81::make_Omega (uint16 nPoint, Eigen::Ref<Eigen::MatrixXd> Omega) {
//was:Mat2<6,9> &B
	for (uint16 i=0; i < 3; i++) {
    Omega(0,i*3+0) = O[nPoint][0+i*3];//exx
    Omega(1,i*3+0) = O[nPoint][1+i*3];//exy
    Omega(1,i*3+1) = O[nPoint][0+i*3];//exy
    Omega(2,i*3+0) = O[nPoint][2+i*3];//exz
    Omega(2,i*3+2) = O[nPoint][0+i*3];//exz
    Omega(3,i*3+1) = O[nPoint][1+i*3];//eyy
    Omega(4,i*3+1) = O[nPoint][2+i*3];//eyz
    Omega(4,i*3+2) = O[nPoint][1+i*3];//eyz
    Omega(5,i*3+2) = O[nPoint][2+i*3];//ezz
	}
}

void ElementSOLID81::getScalar(double& scalar, el_component code, uint16 gp, const double scale) {
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
      error("ElementSOLID81::getScalar: no data for code %d", code);
  }
}


void ElementSOLID81::getVector(double* vector, el_vector code, uint16 gp, const double scale) {
  if (gp == GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 nPoint = 0; nPoint < npow(n_int(),n_dim()); nPoint ++) {
      dWt = g_weight(nPoint);
      getVector(vector, code, nPoint, dWt/dWtSum*scale );
    }
    return;
  }

  Vec<6> vecC;
  double I1C;
  double I2Cz;
  double J;
  switch (code) {
    case VECTOR_IC:
      // be aware about components ordering in C
      vecC = C[gp];
      I1C = vecC[M_XX]+vecC[M_YY]+vecC[M_ZZ]; 
      //I2C with star!
      I2Cz= vecC[M_XX]*vecC[M_YY] + vecC[M_YY]*vecC[M_ZZ] + vecC[M_XX]*vecC[M_ZZ]
                         - vecC[M_XY]*vecC[M_XY] - vecC[M_YZ]*vecC[M_YZ] - vecC[M_XZ]*vecC[M_XZ];
      J = Material::getJ(vecC.ptr());
      vector[0] += I1C*scale;
      vector[1] += I2Cz*scale;
      vector[2] += J*scale;
      break;
    default:
      error("ElementSOLID80::getVector: no data for code %d", code);
  }
}

//return a tensor in a global coordinate system
void  ElementSOLID81::getTensor(MatSym<3>& tensor, el_tensor code, uint16 gp, const double scale) {
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
      error("ElementSOLID81::getTensor: no data for code %d", code);
	}
}
void ElementSOLID81::assemble4(Eigen::Ref<Eigen::MatrixXd> Ke, Eigen::Ref<Eigen::MatrixXd> Fe) {
	assert (Element::n_nodes()*Node::n_dofs() + Element::n_dofs() == Ke.cols());
	assert (Element::n_dofs() == 1);

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
      for (uint16 j=i; j < Element::n_nodes(); j++)
				for (uint16 dj=0; dj < Node::n_dofs(); dj++) {
					if ((i==j) && (dj<di)) {
            continue;
          } else {
						storage->Kij_add(nodes[i],di,nodes[j],dj, Ke(i*3+di,j*3+dj));
					}
				}
	//upper diagonal process for nodes-el dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for(uint16 di=0; di < Node::n_dofs(); di++)
				storage->Kij_add(nodes[i],di, -(int32)getElNum(), 0, Ke(i*3+di, 24));
	//upper diagonal process for el-el dofs
			storage->Kij_add(-(int32)getElNum(), 0, -(int32)getElNum(), 0, Ke(24,24));

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
			storage->Fi_add(nodes[i],di, Fe(i*3+di,0));
  storage->Fi_add(-(int32)getElNum(), 0, Fe(24,0));
}


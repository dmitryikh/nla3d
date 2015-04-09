// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/PLANE41.h"

namespace nla3d {
using namespace math;
using namespace solidmech;

const solidmech::tensorComponents ElementPLANE41::components[3] = {M_XX, M_YY, M_XY};
const uint16 ElementPLANE41::num_components = 3;

//------------------ElementPLANE41--------------------
void ElementPLANE41::pre() {
	S.assign(npow(n_int(),n_dim()), Vec<3>(0.0f, 0.0f, 0.0f));
	C.assign(npow(n_int(),n_dim()), Vec<3>(1.0f, 1.0f, 0.0f));
	O.assign(npow(n_int(),n_dim()), Vec<4>(0.0f, 0.0f, 0.0f, 0.0f));
	if (det.size()==0) {
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->getElementNodes(getElNum(), nodes_p);	
		make_Jacob(getElNum(), nodes_p);
	}

  // register element equations
  for (uint16 i = 0; i < Element::n_nodes(); i++) {
    storage->registerNodeDof(getNodeNumber(i), Dof::UX);
    storage->registerNodeDof(getNodeNumber(i), Dof::UY);
  }
  storage->registerElementDof(getElNum(), Dof::HYDRO_PRESSURE);
}

void ElementPLANE41::build () {
//построение матрицы жесткости и вектора нагрузки для элемента
	Mat<8,8> Kuu; //матрица жесткости элемента
	Mat<8,1> Kup;
	Mat<9,9> Ke;
	double Kpp = 0.0;
	double Fp = 0.0;
	Vec<8> Qe; //вектор узловых сил элемента
	Vec<9> Fe;
  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
	MatSym<3> matD_d;
	Vec<3> vecD_p;
	double p_e = storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE);
  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
  CHECK_NOTNULL(mat);

	double k = mat->getK();
	double dWt; //Gaussian quadrature
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++) {
		dWt = g_weight(nPoint);
    // all meterial functions are waiting [C] for 3D case. So we need to use CVec here.
    CVec[M_XX] = C[nPoint][0];
    CVec[M_YY] = C[nPoint][1];
    CVec[M_XY] = C[nPoint][2];
		mat->getDdDp_UP(num_components, components, CVec.ptr(), p_e, matD_d.ptr(), vecD_p.ptr());
    //matD_d will be 3x3 symmetric matrix. We need to convert it onto 3x3 usual matrix
    Mat<3,3> matE_c = matD_d.toMat();
    Mat<3,1> matE_p;
    matE_p[0][0] = vecD_p[0];
    matE_p[1][0] = vecD_p[1];
    matE_p[2][0] = vecD_p[2];
		double J = solidmech::J_C(CVec.ptr());

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
		Kuu += (matB.transpose() * matE_c * matB * 2.0 + matBomega.transpose() * matS * matBomega)*dWt;
		Fp += (J - 1 - p_e/k)*dWt;
		Qe += (matB.transpose() * S[nPoint] * dWt);
		Kup+= matB.transpose()*matE_p *dWt;
		Kpp -= 1.0/k*dWt;

	}// loop over intergration points
	

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
	assemble(Ke, Fe);
}
//
inline Mat<3,8> ElementPLANE41::make_B(uint16 nPoint) {
	Mat<3,8> B = Mat<3,8>(NjXi[nPoint][0][0], 0.0f, NjXi[nPoint][0][1], 0.0f, NjXi[nPoint][0][2], 0.0f, NjXi[nPoint][0][3], 0.0f,
									0.0f, NjXi[nPoint][1][0], 0.0f, NjXi[nPoint][1][1], 0.0f, NjXi[nPoint][1][2], 0.0f, NjXi[nPoint][1][3],
									NjXi[nPoint][1][0], NjXi[nPoint][0][0], NjXi[nPoint][1][1], NjXi[nPoint][0][1], NjXi[nPoint][1][2], NjXi[nPoint][0][2], NjXi[nPoint][1][3], NjXi[nPoint][0][3]);
	return B;
}
//
Mat<4,8> ElementPLANE41::make_Bomega(uint16 nPoint)
{
	Mat<4,8> Bomega = Mat<4,8>(NjXi[nPoint][0][0], 0.0f, NjXi[nPoint][0][1], 0.0f, NjXi[nPoint][0][2], 0.0f, NjXi[nPoint][0][3], 0.0f,
										NjXi[nPoint][1][0], 0.0f, NjXi[nPoint][1][1], 0.0f, NjXi[nPoint][1][2], 0.0f, NjXi[nPoint][1][3], 0.0f,
										0.0f, NjXi[nPoint][0][0], 0.0f, NjXi[nPoint][0][1], 0.0f, NjXi[nPoint][0][2], 0.0f, NjXi[nPoint][0][3],
										0.0f, NjXi[nPoint][1][0], 0.0f, NjXi[nPoint][1][1], 0.0f, NjXi[nPoint][1][2], 0.0f, NjXi[nPoint][1][3]);
	return Bomega;
}

void ElementPLANE41::update() {
	// get nodal solutions from storage
	Vec<8> U;
  for (uint16 i = 0; i < Element::n_nodes(); i++) {
    U[i*2 + 0] = storage->getDofSolution(getNodeNumber(i), Dof::UX);
    U[i*2 + 1] = storage->getDofSolution(getNodeNumber(i), Dof::UY);
  }
  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
  CHECK_NOTNULL(mat);
  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
	//восстанавливаем преращение давления
	double p_e = storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE);

  for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++) {
		Mat<4,8> matBomega = make_Bomega(nPoint);
		O[nPoint] = matBomega * U;
		C[nPoint][0] = 1.0f + 2*O[nPoint][0]+1.0f*(O[nPoint][0]*O[nPoint][0]+O[nPoint][2]*O[nPoint][2]);
		C[nPoint][1]=1.0f + 2*O[nPoint][3]+1.0f*(O[nPoint][3]*O[nPoint][3]+O[nPoint][1]*O[nPoint][1]);
		C[nPoint][2]=O[nPoint][1]+O[nPoint][2]+O[nPoint][0]*O[nPoint][1]+O[nPoint][2]*O[nPoint][3];
		//восстановление напряжений Пиолы-Кирхгоффа из текущего состояния
    //all meterial functions are waiting [C] for 3D case. So we need to use CVec here.
    CVec[M_XX] = C[nPoint][0];
    CVec[M_YY] = C[nPoint][1];
    CVec[M_XY] = C[nPoint][2];
		mat->getS_UP (num_components, components, CVec.ptr(), p_e, S[nPoint].ptr());
	}
}

void ElementPLANE41::getScalar(double& scalar, query::scalarQuery code, uint16 gp, const double scale) {
	//see codes in sys.h
	//gp - needed gauss point 
  if (gp == query::GP_MEAN) { //need to average result over the element
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
    case query::SCALAR_SP:
      scalar += storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE) * scale;
			break;
    default:
      LOG_N_TIMES(10, WARNING) << "No data for code " << code;
  }
}

void ElementPLANE41::getVector(double* vector, query::vectorQuery code, uint16 gp, const double scale) {
  if (gp == query::GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 nPoint = 0; nPoint < npow(n_int(),n_dim()); nPoint ++) {
      dWt = g_weight(nPoint);
      getVector(vector, code, nPoint, dWt/dWtSum*scale );
    }
    return;
  }
  double IC[3];
  Vec<6> CVec;
  switch (code) {
    case query::VECTOR_IC:
      CVec[M_XZ] = 0.0;
      CVec[M_YZ] = 0.0;
      CVec[M_ZZ] = 1.0;
      CVec[M_XX] = C[gp][0];
      CVec[M_YY] = C[gp][1];
      CVec[M_XY] = C[gp][2];
      solidmech::IC_C(CVec.ptr(), IC);
      vector[0] += IC[0] * scale;
      vector[1] += IC[1] * scale;
      vector[2] += IC[2] * scale;
      break;
    default:
      LOG_N_TIMES(10, WARNING) << "No data for code " << code;
  }
}
//return a tensor in a global coordinate system
void  ElementPLANE41::getTensor(math::MatSym<3>& tensor, query::tensorQuery code, uint16 gp, const double scale) {
  if (gp == query::GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 nPoint = 0; nPoint < npow(n_int(),n_dim()); nPoint ++) {
      dWt = g_weight(nPoint);
      getTensor(tensor, code, nPoint, dWt/dWtSum*scale);
    }
    return;
  }
	assert (npow(n_int(),n_dim()));

  MatSym<3> matS;
  Mat_Hyper_Isotrop_General* mat;
  Mat2<3,3> matF;
  double J;

  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
  CVec[M_XX] = C[gp][0];
  CVec[M_YY] = C[gp][1];
  CVec[M_XY] = C[gp][2];

  double p_e;

	switch (code) {
    case query::TENSOR_COUCHY:

      matF.zeros();
      matF.data[0][0] = 1+O[gp][0]; //11
      matF.data[0][1] = O[gp][1];  //12
      matF.data[1][0] = O[gp][2];  //21
      matF.data[1][1] = 1+O[gp][3];//22
      matF.data[2][2] = 1; //33

      J = matF.data[0][0]*(matF.data[1][1]*matF.data[2][2]-matF.data[1][2]*matF.data[2][1])-matF.data[0][1]*(matF.data[1][0]*matF.data[2][2]-matF.data[1][2]*matF.data[2][0])+matF.data[0][2]*(matF.data[1][0]*matF.data[2][1]-matF.data[1][1]*matF.data[2][0]);
      //In order to complete matS (3x3 symmetric matrix, PK2 tensor) we need 
      //to know S33 component: 
      //1) One solution is to calculate S33 on every solution step
      //and store it in S[nPoint] vector.
      //2) Second solution is to resotre S33 right here.
      //Now 2) is working.
      p_e = storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE);
      mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
      CHECK_NOTNULL(mat);
      mat->getS_UP (6, solidmech::defaultTensorComponents, CVec.ptr(), p_e, matS.data);
      matBTDBprod (matF, matS, 1.0/J, tensor); //Symmetric Couchy tensor
      break;
    case query::TENSOR_PK2:
      mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
      CHECK_NOTNULL(mat);
      p_e = storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE);
      mat->getS_UP (6, solidmech::defaultTensorComponents, CVec.ptr(), p_e, matS.data);
      tensor.data[0] += matS.data[0]*scale;
      tensor.data[1] += matS.data[1]*scale;
      tensor.data[2] += matS.data[2]*scale;
      tensor.data[3] += matS.data[3]*scale;
      tensor.data[4] += matS.data[4]*scale;
      tensor.data[5] += matS.data[5]*scale;
      break;
    case query::TENSOR_C:
      tensor.data[0] += CVec[M_XX]*scale;
      tensor.data[1] += CVec[M_XY]*scale;
      tensor.data[2] += CVec[M_XZ]*scale;
      tensor.data[3] += CVec[M_YY]*scale;
      tensor.data[4] += CVec[M_YZ]*scale;
      tensor.data[5] += CVec[M_ZZ]*scale;
      break;
    case query::TENSOR_E:
      tensor.data[0] += (CVec[M_XX]-1.0)*0.5*scale;
      tensor.data[1] += CVec[M_XY]*0.5*scale;
      tensor.data[2] += CVec[M_XZ]*0.5*scale;
      tensor.data[3] += (CVec[M_YY]-1.0)*0.5*scale;
      tensor.data[4] += CVec[M_YZ]*0.5*scale;
      tensor.data[5] += (CVec[M_ZZ]-1.0)*0.5*scale;
      break;
    default:
      LOG_N_TIMES(10, WARNING) << "No data for code " << code;
	}
}

} // namespace nla3d

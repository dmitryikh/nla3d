// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/SOLID81.h"

namespace nla3d {
using namespace math;
using namespace solidmech;

void ElementSOLID81::pre() {
	S.assign(npow(n_int(),n_dim()), Vec<6>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
	C.assign(npow(n_int(),n_dim()), Vec<6>(1.0, 0.0, 0.0, 1.0, 0.0, 1.0));
	O.assign(npow(n_int(),n_dim()), Vec<9>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

	if (det.size()==0) {
		Node** nodes_p = new Node*[Element::n_nodes()];
		storage->getElementNodes(getElNum(), nodes_p);	
		make_Jacob(getElNum(), nodes_p);
    delete[] nodes_p;
	}

  // register element equations
  for (uint16 i = 0; i < Element::n_nodes(); i++) {
    storage->registerNodeDof(getNodeNumber(i), Dof::UX);
    storage->registerNodeDof(getNodeNumber(i), Dof::UY);
    storage->registerNodeDof(getNodeNumber(i), Dof::UZ);
  }
  storage->registerElementDof(getElNum(), Dof::HYDRO_PRESSURE);
}


void ElementSOLID81::build()
{
//построение матрицы жесткости и вектора нагрузки для элемента
	double Kpp = 0.0;
	double Fp = 0.0;
	
	Vec<24> Kup;
	Vec<24> Fu; //вектор узловых сил элемента
	Vec<24> F_ext; //вектор внешних сил (пока не подсчитывается)
  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
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
	double p_e = storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE);
	double dWt; //множитель при суммировании квадратур Гаусса
	Kuu.zeros();
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++) {
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		dWt = g_weight(nPoint);

		mat->getDdDp_UP(6, solidmech::defaultTensorComponents, C[nPoint].ptr(), p_e, matD_d.ptr(), vecD_p.ptr());
		double J = solidmech::J_C(C[nPoint].ptr());
		matB.zeros();
		matS.zeros();
		matO.zeros();
		matB_NL.zeros();

		make_B_L(nPoint, matB);
		make_S(nPoint, matS); //matrix 9x9 with 3x3 stress tenros blocks
		make_Omega(nPoint, matO);
		make_B_NL(nPoint, matB_NL);
    // matB = matB + matO * matB_NL * 2
		matABprod(matO, matB_NL, 2.0, matB);
    // Kuu = Kuu + (matB^T * matD_d * matB) * 0.5*dWt
		matBTDBprod(matB, matD_d, 0.5*dWt, Kuu);
    // Kuu = Kuu + (matB_NL^T * matS * matB_NL) * dWt
		matBTDBprod(matB_NL, matS, dWt, Kuu);

		// Fu = Fu +  matB^T * S[nPoint] * (-0.5*dWt);
		matBTVprod(matB,S[nPoint], -0.5*dWt, Fu);

		// Kup = Kup +  matB^T * vecD_p * (dWt*0.5);
		matBTVprod(matB,vecD_p, 0.5*dWt, Kup);

		Fp += -(J - 1 - p_e/k)*dWt;
		Kpp += -1.0/k*dWt;
	}//прошлись по всем точкам интегрирования
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);
	assemble3(Kuu, Kup, Kpp, Fu,Fp);
}


void ElementSOLID81::update()
{
	// get nodal solutions from storage
	Vec<24> U;
  for (uint16 i = 0; i < Element::n_nodes(); i++) {
    U[i*3 + 0] = storage->getDofSolution(getNodeNumber(i), Dof::UX);
    U[i*3 + 1] = storage->getDofSolution(getNodeNumber(i), Dof::UY);
    U[i*3 + 2] = storage->getDofSolution(getNodeNumber(i), Dof::UZ);
  }
	Mat2<9,24> B_NL;
	Vec<6> vecC;
	double p_e = storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE);

  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
  if (mat == NULL) {
    error("ElementSOLID81::update: material is not derived from Mat_Hyper_Isotrop_General");
  }
	for (uint16 nPoint=0; (int32) nPoint < npow(Element::n_int(),n_dim()); nPoint++)
	{
		GlobStates.setuint16(States::UI16_CURINTPOINT, nPoint);
		B_NL.zeros();
		make_B_NL(nPoint, B_NL);
		O[nPoint].zeros();
    // O[nPoint] = B_NL * U
		matBVprod(B_NL, U, 1.0, O[nPoint]);
		C[nPoint][M_XX] = 1.0+2*O[nPoint][0]+pow(O[nPoint][0],2)+pow(O[nPoint][3],2)+pow(O[nPoint][6],2);	//C11
		C[nPoint][M_YY] = 1.0+2*O[nPoint][4]+pow(O[nPoint][1],2)+pow(O[nPoint][4],2)+pow(O[nPoint][7],2);	//C22
		C[nPoint][M_ZZ] = 1.0+2*O[nPoint][8]+pow(O[nPoint][2],2)+pow(O[nPoint][5],2)+pow(O[nPoint][8],2);	//C33
		C[nPoint][M_XY] = O[nPoint][1]+O[nPoint][3]+O[nPoint][0]*O[nPoint][1]+O[nPoint][3]*O[nPoint][4]+O[nPoint][6]*O[nPoint][7];  //C12
		C[nPoint][M_YZ] = O[nPoint][5]+O[nPoint][7]+O[nPoint][1]*O[nPoint][2]+O[nPoint][4]*O[nPoint][5]+O[nPoint][7]*O[nPoint][8];	//C23
		C[nPoint][M_XZ] = O[nPoint][2]+O[nPoint][6]+O[nPoint][0]*O[nPoint][2]+O[nPoint][3]*O[nPoint][5]+O[nPoint][6]*O[nPoint][8];	//C13

		// get new PK2 stresses from update C tensor
		mat->getS_UP (6, solidmech::defaultTensorComponents, C[nPoint].ptr(), p_e, S[nPoint].ptr());
	}
	GlobStates.undefineuint16(States::UI16_CURINTPOINT);
	GlobStates.undefinedouble(States::DOUBLE_HYDPRES);
}


void ElementSOLID81::make_B_L (uint16 nPoint, Mat2<6,24> &B)
{
	double *B_L = B.ptr();
	for (uint16 i=0; i < 8; i++) {
		B_L[0*24+(i*3+0)] += 2*NjXi[nPoint][0][i];//exx
		B_L[1*24+(i*3+0)] += 2*NjXi[nPoint][1][i];//exy
		B_L[1*24+(i*3+1)] += 2*NjXi[nPoint][0][i];//exy
		B_L[2*24+(i*3+0)] += 2*NjXi[nPoint][2][i];//exz
		B_L[2*24+(i*3+2)] += 2*NjXi[nPoint][0][i];//exz
		B_L[3*24+(i*3+1)] += 2*NjXi[nPoint][1][i];//eyy
		B_L[4*24+(i*3+1)] += 2*NjXi[nPoint][2][i];//eyz
		B_L[4*24+(i*3+2)] += 2*NjXi[nPoint][1][i];//eyz
		B_L[5*24+(i*3+2)] += 2*NjXi[nPoint][2][i];//ezz
	}
}


void ElementSOLID81::make_B_NL (uint16 nPoint,  Mat2<9,24> &B)
{
	double *B_NL = B.ptr();
	for (uint16 i=0; i < 8; i++) {
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



void ElementSOLID81::make_S (uint16 nPoint, MatSym<9> &B)
{
	double *Sp = B.ptr();
	Sp[0]	+= S[nPoint][M_XX];
	Sp[1]	+= S[nPoint][M_XY];
	Sp[2]	+= S[nPoint][M_XZ];
	Sp[9]	+= S[nPoint][M_YY];
	Sp[10]+= S[nPoint][M_YZ];
	Sp[17]+= S[nPoint][M_ZZ];

	Sp[24]+= S[nPoint][M_XX];
	Sp[25]+= S[nPoint][M_XY];
	Sp[26]+= S[nPoint][M_XZ];
	Sp[30]+= S[nPoint][M_YY];
	Sp[31]+= S[nPoint][M_YZ];
	Sp[35]+= S[nPoint][M_ZZ];

	Sp[39]	+= S[nPoint][M_XX];
	Sp[40]	+= S[nPoint][M_XY];
	Sp[41]	+= S[nPoint][M_XZ];
	Sp[42]	+= S[nPoint][M_YY];
	Sp[43]	+= S[nPoint][M_YZ];
  Sp[44]  += S[nPoint][M_ZZ];
}


void ElementSOLID81::make_Omega (uint16 nPoint, Mat2<6,9> &B)
{
	double *Omega = B.ptr();
	for (uint16 i=0; i < 3; i++) {
    Omega[0*9+(i*3+0)] = O[nPoint][0+i*3];//exx
    Omega[1*9+(i*3+0)] = O[nPoint][1+i*3];//exy
    Omega[1*9+(i*3+1)] = O[nPoint][0+i*3];//exy
    Omega[2*9+(i*3+0)] = O[nPoint][2+i*3];//exz
    Omega[2*9+(i*3+2)] = O[nPoint][0+i*3];//exz
    Omega[3*9+(i*3+1)] = O[nPoint][1+i*3];//eyy
    Omega[4*9+(i*3+1)] = O[nPoint][2+i*3];//eyz
    Omega[4*9+(i*3+2)] = O[nPoint][1+i*3];//eyz
    Omega[5*9+(i*3+2)] = O[nPoint][2+i*3];//ezz
	}
}

void ElementSOLID81::getScalar(double& scalar, query::scalarQuery code, uint16 gp, const double scale) {
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
	assert (gp < npow(n_int(),n_dim()));
  double tmp[3];
  Mat_Hyper_Isotrop_General* mat;
  double J;
  switch (code) {
    case query::SCALAR_SP:
      scalar += storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE) * scale;
			break;
    case query::SCALAR_WU:
      tmp[0] = 0.0;
      tmp[1] = 0.0;
      tmp[2] = 0.0;
      getVector(tmp, query::VECTOR_IC, gp, 1.0);
      mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
      if (mat == NULL) {
        error("ElementSOLID81::getScalar: material is not derived from Mat_Hyper_Isotrop_General");
      }
      scalar += mat->W(tmp[0], tmp[1], tmp[2])*scale;
      break;
    case query::SCALAR_WP:
      J = solidmech::J_C(C[gp].ptr());
      mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
      if (mat == NULL) {
        error("ElementSOLID81::getScalar: material is not derived from Mat_Hyper_Isotrop_General");
      }
      scalar += 0.5*mat->getK()*(J-1.0)*(J-1.0)*scale;
      break;
    case query::SCALAR_VOL:
      scalar += volume()*scale;
      break;
    default:
      error("ElementSOLID81::getScalar: no data for code %d", code);
  }
}


void ElementSOLID81::getVector(double* vector, query::vectorQuery code, uint16 gp, const double scale) {
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
  switch (code) {
    case query::VECTOR_IC:
      solidmech::IC_C(C[gp].ptr(), IC);
      vector[0] += IC[0] * scale;
      vector[1] += IC[1] * scale;
      vector[2] += IC[2] * scale;
      break;
    default:
      error("ElementSOLID81::getVector: no data for code %d", code);
  }
}



//return a tensor in a global coordinate system
void  ElementSOLID81::getTensor(math::MatSym<3>& tensor, query::tensorQuery code, uint16 gp, const double scale) {
  if (gp == query::GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 nPoint = 0; nPoint < npow(n_int(),n_dim()); nPoint ++) {
      dWt = g_weight(nPoint);
      getTensor(tensor, code, nPoint, dWt/dWtSum*scale);
    }
    return;
  }
	assert (gp < npow(n_int(),n_dim()));

  Mat2<3,3> matF;
  MatSym<3> matS;
  double J;
  double cInv[6];
  double pe;

	switch (code) {
    case query::TENSOR_COUCHY:

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

      J = solidmech::J_C(C[gp].ptr());

      //deviatoric part of S: Sd = S[gp]
      //hydrostatic part of S: Sp = p * J * C^(-1)
      solidmech::invC_C (C[gp].ptr(), J, cInv); 
	    pe = storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE);
      // TODO: it seems that S[gp] contains deviatoric + pressure already..
      // we dont need to sum it again
      for (size_t i = 0; i < 6; i++) {
        matS.data[i] = S[gp][i] + pe * J * cInv[i];
      }
      matBTDBprod (matF, matS, 1.0/J*scale, tensor); //Symmetric Couchy tensor
      break;
    case query::TENSOR_PK2:
      // hydrostatic part of S: Sp = p * J * C^(-1)
      J = solidmech::J_C(C[gp].ptr());
      solidmech::invC_C (C[gp].ptr(), J, cInv); 
	    pe = storage->getDofSolution(- static_cast<int32> (getElNum()), Dof::HYDRO_PRESSURE);
      for (size_t i = 0; i < 6; i++) {
        tensor.data[i] += (S[gp][i] + pe * J * cInv[i]) * scale;
      }
      break;
    case query::TENSOR_C:
      tensor.data[0] += C[gp][M_XX]*scale;
      tensor.data[1] += C[gp][M_XY]*scale;
      tensor.data[2] += C[gp][M_XZ]*scale;
      tensor.data[3] += C[gp][M_YY]*scale;
      tensor.data[4] += C[gp][M_YZ]*scale;
      tensor.data[5] += C[gp][M_ZZ]*scale;
      break;
    case query::TENSOR_E:
      tensor.data[0] += (C[gp][M_XX]-1.0)*0.5*scale;
      tensor.data[1] += C[gp][M_XY]*0.5*scale;
      tensor.data[2] += C[gp][M_XZ]*0.5*scale;
      tensor.data[3] += (C[gp][M_YY]-1.0)*0.5*scale;
      tensor.data[4] += C[gp][M_YZ]*0.5*scale;
      tensor.data[5] += (C[gp][M_ZZ]-1.0)*0.5*scale;
      break;
    default:
      error("ElementSOLID81::getTensor: no data for code %d", code);
	}
}

} // namespace nla3d

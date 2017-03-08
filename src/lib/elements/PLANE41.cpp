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
  if (det.size()==0) {
    makeJacob();
  }

  S.assign(nOfIntPoints(), Vec<3>(0.0f, 0.0f, 0.0f));
  C.assign(nOfIntPoints(), Vec<3>(1.0f, 1.0f, 0.0f));
  O.assign(nOfIntPoints(), Vec<4>(0.0f, 0.0f, 0.0f, 0.0f));

  // register element equations
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY});
  }
  storage->addElementDof(getElNum(), {Dof::HYDRO_PRESSURE});
}

void ElementPLANE41::buildK() {
  Mat<8,8> Kuu;  // displacement stiff. matrix
  Mat<8,1> Kup;
  Mat<9,9> Ke;  // element stiff. matrix
  double Kpp = 0.0;
  double Fp = 0.0;
  Vec<8> Qe;  // displacement rhs
  Vec<9> Fe;  // element
  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
  MatSym<3> matD_d;
  Vec<3> vecD_p;
  double p_e = storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE);
  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
  CHECK_NOTNULL(mat);

  double k = mat->getK();
  double dWt; //Gaussian quadrature
  for (uint16 np=0; np < nOfIntPoints(); np++) {
    dWt = intWeight(np);
    // all meterial functions are waiting [C] for 3D case. So we need to use CVec here.
    CVec[M_XX] = C[np][0];
    CVec[M_YY] = C[np][1];
    CVec[M_XY] = C[np][2];
    mat->getDdDp_UP(num_components, components, CVec.ptr(), p_e, matD_d.ptr(), vecD_p.ptr());
    //matD_d will be 3x3 symmetric matrix. We need to convert it onto 3x3 usual matrix
    Mat<3,3> matE_c = matD_d.toMat();
    Mat<3,1> matE_p;
    matE_p[0][0] = vecD_p[0];
    matE_p[1][0] = vecD_p[1];
    matE_p[2][0] = vecD_p[2];
    double J = solidmech::J_C(CVec.ptr());

    Mat<3,8> matB = make_B(np);
    //матрица S для матричного умножения
    Mat<4,4> matS = Mat<4,4>(S[np][0],S[np][2], 0.0, 0.0, 
                S[np][2],S[np][1], 0.0, 0.0,
                0.0, 0.0, S[np][0],S[np][2],
                0.0, 0.0, S[np][2],S[np][1]);
    //матрица Омега.используется для составления 
    //матр. накопленных линейных деформаций к текущему шагу
    Mat<3,4> matO = Mat<3,4>(O[np][0], 0.0, O[np][2], 0.0,
                0.0, O[np][1], 0.0, O[np][3],
                O[np][1], O[np][0], O[np][3], O[np][2]);

    Mat<4,8> matBomega = make_Bomega(np);
    Mat<3,8> matBl = matO * matBomega;
    matB += matBl;
    Kuu += (matB.transpose() * matE_c * matB * 2.0 + matBomega.transpose() * matS * matBomega)*dWt;
    Fp += (J - 1 - p_e/k)*dWt;
    Qe += (matB.transpose() * S[np] * dWt);
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
  assembleK(Ke, Fe);
}
//
inline Mat<3,8> ElementPLANE41::make_B(uint16 np) {
  Mat<3,8> B = Mat<3,8>(NiXj[np][0][0], 0.0f, NiXj[np][1][0], 0.0f, NiXj[np][2][0], 0.0f, NiXj[np][3][0], 0.0f,
                  0.0f, NiXj[np][0][1], 0.0f, NiXj[np][1][1], 0.0f, NiXj[np][2][1], 0.0f, NiXj[np][3][1],
                  NiXj[np][0][1], NiXj[np][0][0], NiXj[np][1][1], NiXj[np][1][0], NiXj[np][2][1], NiXj[np][2][0], NiXj[np][3][1], NiXj[np][3][0]);
  return B;
}
//
Mat<4,8> ElementPLANE41::make_Bomega(uint16 np) {
  Mat<4,8> Bomega = Mat<4,8>(NiXj[np][0][0], 0.0f, NiXj[np][1][0], 0.0f, NiXj[np][2][0], 0.0f, NiXj[np][3][0], 0.0f,
                    NiXj[np][0][1], 0.0f, NiXj[np][1][1], 0.0f, NiXj[np][2][1], 0.0f, NiXj[np][3][1], 0.0f,
                    0.0f, NiXj[np][0][0], 0.0f, NiXj[np][1][0], 0.0f, NiXj[np][2][0], 0.0f, NiXj[np][3][0],
                    0.0f, NiXj[np][0][1], 0.0f, NiXj[np][1][1], 0.0f, NiXj[np][2][1], 0.0f, NiXj[np][3][1]);
  return Bomega;
}

void ElementPLANE41::update() {
  // get nodal solutions from storage
  Vec<8> U;
  for (uint16 i = 0; i < getNNodes(); i++) {
    U[i*2 + 0] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U[i*2 + 1] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
  }
  Mat_Hyper_Isotrop_General* mat = dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial());
  CHECK_NOTNULL(mat);
  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
  //восстанавливаем преращение давления
  double p_e = storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE);

  for (uint16 np=0; np < nOfIntPoints(); np++) {
    Mat<4,8> matBomega = make_Bomega(np);
    O[np] = matBomega * U;
    C[np][0] = 1.0f + 2*O[np][0]+1.0f*(O[np][0]*O[np][0]+O[np][2]*O[np][2]);
    C[np][1]=1.0f + 2*O[np][3]+1.0f*(O[np][3]*O[np][3]+O[np][1]*O[np][1]);
    C[np][2]=O[np][1]+O[np][2]+O[np][0]*O[np][1]+O[np][2]*O[np][3];
    //восстановление напряжений Пиолы-Кирхгоффа из текущего состояния
    //all meterial functions are waiting [C] for 3D case. So we need to use CVec here.
    CVec[M_XX] = C[np][0];
    CVec[M_YY] = C[np][1];
    CVec[M_XY] = C[np][2];
    mat->getS_UP (num_components, components, CVec.ptr(), p_e, S[np].ptr());
  }
}


bool ElementPLANE41::getScalar(double* scalar, scalarQuery query, uint16 gp, const double scale) {
  //see queries in query.h
  //gp - needed gauss point 
  assert(scalar != nullptr);

  if (gp == GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 np = 0; np < nOfIntPoints(); np ++) {
      dWt = intWeight(np);
      bool ret = getScalar(scalar, query, np, dWt / dWtSum * scale);
      if(ret == false) return false;
    }
    return true;
  }

  switch (query) {
    case scalarQuery::SP:
      *scalar += storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE) * scale;
      return true;
  }
  return false;
}


bool ElementPLANE41::getVector(math::Vec<3>* vector, vectorQuery query, uint16 gp, const double scale) {
  assert(vector != nullptr);

  if (gp == GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 np = 0; np < nOfIntPoints(); np ++) {
      dWt = intWeight(np);
      bool ret = getVector(vector, query, np, dWt / dWtSum * scale);
      if(ret == false) return false;
    }
    return true;
  }

  double IC[3];
  Vec<6> CVec;
  switch (query) {
    case vectorQuery::IC:
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
      return true;
  }
  return false;
}


//return a tensor in a global coordinate system
bool  ElementPLANE41::getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale) {
  assert(tensor != nullptr);
  if (gp == GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 np = 0; np < nOfIntPoints(); np ++) {
      dWt = intWeight(np);
      bool ret = getTensor(tensor, query, np, dWt / dWtSum * scale);
      if(ret == false) return false;
    }
    return true;
  }

  // here we obtain result for query on particular Gaussian point
  assert (gp < nOfIntPoints());

  MatSym<3> matS;
  Mat_Hyper_Isotrop_General* mat;
  Mat<3,3> matF;
  double J;

  Vec<6> CVec;
  CVec[M_XZ] = 0.0;
  CVec[M_YZ] = 0.0;
  CVec[M_ZZ] = 1.0;
  CVec[M_XX] = C[gp][0];
  CVec[M_YY] = C[gp][1];
  CVec[M_XY] = C[gp][2];

  double p_e;

  switch (query) {
    case tensorQuery::COUCHY:

      matF.zero();
      matF.data[0][0] = 1+O[gp][0]; //11
      matF.data[0][1] = O[gp][1];  //12
      matF.data[1][0] = O[gp][2];  //21
      matF.data[1][1] = 1+O[gp][3];//22
      matF.data[2][2] = 1; //33

      J = matF.data[0][0]*(matF.data[1][1]*matF.data[2][2]-matF.data[1][2]*matF.data[2][1])-
          matF.data[0][1]*(matF.data[1][0]*matF.data[2][2]-matF.data[1][2]*matF.data[2][0])+
          matF.data[0][2]*(matF.data[1][0]*matF.data[2][1]-matF.data[1][1]*matF.data[2][0]);
      //In order to complete matS (3x3 symmetric matrix, PK2 tensor) we need 
      //to know S33 component: 
      //1) One solution is to calculate S33 on every solution step
      //and store it in S[np] vector.
      //2) Second solution is to resotre S33 right here.
      //Now 2) is working.
      p_e = storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE);
      mat = dynamic_cast<Mat_Hyper_Isotrop_General*>(storage->getMaterial());
      CHECK_NOTNULL(mat);
      mat->getS_UP(6, solidmech::defaultTensorComponents, CVec.ptr(), p_e, matS.data);
      matBTDBprod(matF, matS, 1.0/J, *tensor); //Symmetric Couchy tensor
      return true;

    case tensorQuery::PK2:
      mat = dynamic_cast<Mat_Hyper_Isotrop_General*>(storage->getMaterial());
      CHECK_NOTNULL(mat);
      p_e = storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE);
      mat->getS_UP(6, solidmech::defaultTensorComponents, CVec.ptr(), p_e, matS.data);
      tensor->data[0] += matS.data[0]*scale;
      tensor->data[1] += matS.data[1]*scale;
      tensor->data[2] += matS.data[2]*scale;
      tensor->data[3] += matS.data[3]*scale;
      tensor->data[4] += matS.data[4]*scale;
      tensor->data[5] += matS.data[5]*scale;
      return true;

    case tensorQuery::C:
      tensor->data[0] += CVec[M_XX]*scale;
      tensor->data[1] += CVec[M_XY]*scale;
      tensor->data[2] += CVec[M_XZ]*scale;
      tensor->data[3] += CVec[M_YY]*scale;
      tensor->data[4] += CVec[M_YZ]*scale;
      tensor->data[5] += CVec[M_ZZ]*scale;
      return true;

    case tensorQuery::E:
      tensor->data[0] += (CVec[M_XX]-1.0)*0.5*scale;
      tensor->data[1] += CVec[M_XY]*0.5*scale;
      tensor->data[2] += CVec[M_XZ]*0.5*scale;
      tensor->data[3] += (CVec[M_YY]-1.0)*0.5*scale;
      tensor->data[4] += CVec[M_YZ]*0.5*scale;
      tensor->data[5] += (CVec[M_ZZ]-1.0)*0.5*scale;
      return true;
  }
  return false;
}

} // namespace nla3d

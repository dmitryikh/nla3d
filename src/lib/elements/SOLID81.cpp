// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/SOLID81.h"

namespace nla3d {
using namespace math;
using namespace solidmech;


void ElementSOLID81::pre() {
  if (det.size()==0) {
    makeJacob();
  }

  S.assign(nOfIntPoints(), Vec<6>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
  C.assign(nOfIntPoints(), Vec<6>(1.0, 0.0, 0.0, 1.0, 0.0, 1.0));
  O.assign(nOfIntPoints(), Vec<9>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

  // register element equations
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY, Dof::UZ});
  }
  storage->addElementDof(getElNum(), {Dof::HYDRO_PRESSURE});
}


void ElementSOLID81::buildK() {

  double Kpp = 0.0;
  double Fp = 0.0;

  Vec<24> Kup;
  Vec<24> Fu; //вектор узловых сил элемента
  Vec<24> F_ext; //вектор внешних сил (пока не подсчитывается)
  Mat_Hyper_Isotrop_General* mat = CHECK_NOTNULL( dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial()));
  double k = mat->getK();
  MatSym<6> matD_d;
  Vec<6> vecD_p;
  Vec<6> vecC;

  Mat<6,24> matB;
  MatSym<9> matS;
  Mat<6,9> matO;
  Mat<9,24> matB_NL;
  MatSym<24> Kuu; //матрица жесткости перед вектором перемещений
  double p_e = storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE);
  double dWt; //множитель при суммировании квадратур Гаусса
  Kuu.zero();
  for (uint16 np = 0; np < nOfIntPoints(); np++) {
    dWt = intWeight(np);

    mat->getDdDp_UP(6, solidmech::defaultTensorComponents, C[np].ptr(), p_e, matD_d.ptr(), vecD_p.ptr());
    double J = solidmech::J_C(C[np].ptr());
    matB.zero();
    matS.zero();
    matO.zero();
    matB_NL.zero();

    make_B_L(np, matB);
    make_S(np, matS); //matrix 9x9 with 3x3 stress tenros blocks
    make_Omega(np, matO);
    make_B_NL(np, matB_NL);
    // matB = matB + matO * matB_NL * 2
    matABprod(matO, matB_NL, 2.0, matB);
    // Kuu = Kuu + (matB^T * matD_d * matB) * 0.5*dWt
    matBTDBprod(matB, matD_d, 0.5*dWt, Kuu);
    // Kuu = Kuu + (matB_NL^T * matS * matB_NL) * dWt
    matBTDBprod(matB_NL, matS, dWt, Kuu);

    // Fu = Fu +  matB^T * S[np] * (-0.5*dWt);
    matBTVprod(matB,S[np], -0.5*dWt, Fu);

    // Kup = Kup +  matB^T * vecD_p * (dWt*0.5);
    matBTVprod(matB,vecD_p, 0.5*dWt, Kup);

    Fp += -(J - 1 - p_e/k)*dWt;
    Kpp += -1.0/k*dWt;
  }//прошлись по всем точкам интегрирования
  assemble3(Kuu, Kup, Kpp, Fu,Fp);
}


void ElementSOLID81::update()
{
  // get nodal solutions from storage
  Vec<24> U;
  for (uint16 i = 0; i < getNNodes(); i++) {
    U[i*3 + 0] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U[i*3 + 1] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
    U[i*3 + 2] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UZ);
  }
  Mat<9,24> B_NL;
  Vec<6> vecC;
  double p_e = storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE);

  Mat_Hyper_Isotrop_General* mat = CHECK_NOTNULL(dynamic_cast<Mat_Hyper_Isotrop_General*> (storage->getMaterial()));
  for (uint16 np = 0; np < nOfIntPoints(); np++) {
    B_NL.zero();
    make_B_NL(np, B_NL);
    O[np].zero();
    // O[np] = B_NL * U
    matBVprod(B_NL, U, 1.0, O[np]);
    C[np][M_XX] = 1.0+2*O[np][0]+pow(O[np][0],2)+pow(O[np][3],2)+pow(O[np][6],2); //C11
    C[np][M_YY] = 1.0+2*O[np][4]+pow(O[np][1],2)+pow(O[np][4],2)+pow(O[np][7],2); //C22
    C[np][M_ZZ] = 1.0+2*O[np][8]+pow(O[np][2],2)+pow(O[np][5],2)+pow(O[np][8],2); //C33
    C[np][M_XY] = O[np][1]+O[np][3]+O[np][0]*O[np][1]+O[np][3]*O[np][4]+O[np][6]*O[np][7];  //C12
    C[np][M_YZ] = O[np][5]+O[np][7]+O[np][1]*O[np][2]+O[np][4]*O[np][5]+O[np][7]*O[np][8];  //C23
    C[np][M_XZ] = O[np][2]+O[np][6]+O[np][0]*O[np][2]+O[np][3]*O[np][5]+O[np][6]*O[np][8];  //C13

    // get new PK2 stresses from update C tensor
    mat->getS_UP (6, solidmech::defaultTensorComponents, C[np].ptr(), p_e, S[np].ptr());
  }
}


void ElementSOLID81::make_B_L (uint16 np, Mat<6,24> &B)
{
  double *B_L = B.ptr();
  for (uint16 i=0; i < 8; i++) {
    B_L[0*24+(i*3+0)] += 2*NiXj[np][i][0];  // exx
    B_L[1*24+(i*3+0)] += 2*NiXj[np][i][1];  // exy
    B_L[1*24+(i*3+1)] += 2*NiXj[np][i][0];  // exy
    B_L[2*24+(i*3+0)] += 2*NiXj[np][i][2];  // exz
    B_L[2*24+(i*3+2)] += 2*NiXj[np][i][0];  // exz
    B_L[3*24+(i*3+1)] += 2*NiXj[np][i][1];  // eyy
    B_L[4*24+(i*3+1)] += 2*NiXj[np][i][2];  // eyz
    B_L[4*24+(i*3+2)] += 2*NiXj[np][i][1];  // eyz
    B_L[5*24+(i*3+2)] += 2*NiXj[np][i][2];  // ezz
  }
}


void ElementSOLID81::make_B_NL (uint16 np,  Mat<9,24> &B)
{
  double *B_NL = B.ptr();
  for (uint16 i=0; i < 8; i++) {
    B_NL[0*24+(i*3+0)] += NiXj[np][i][0];
    B_NL[1*24+(i*3+0)] += NiXj[np][i][1];
    B_NL[2*24+(i*3+0)] += NiXj[np][i][2];
    B_NL[3*24+(i*3+1)] += NiXj[np][i][0];
    B_NL[4*24+(i*3+1)] += NiXj[np][i][1];
    B_NL[5*24+(i*3+1)] += NiXj[np][i][2];
    B_NL[6*24+(i*3+2)] += NiXj[np][i][0];
    B_NL[7*24+(i*3+2)] += NiXj[np][i][1];
    B_NL[8*24+(i*3+2)] += NiXj[np][i][2];
  }
}



void ElementSOLID81::make_S (uint16 np, MatSym<9> &B) {
  double *Sp = B.ptr();
  Sp[0]  += S[np][M_XX];
  Sp[1]  += S[np][M_XY];
  Sp[2]  += S[np][M_XZ];
  Sp[9]  += S[np][M_YY];
  Sp[10] += S[np][M_YZ];
  Sp[17] += S[np][M_ZZ];

  Sp[24] += S[np][M_XX];
  Sp[25] += S[np][M_XY];
  Sp[26] += S[np][M_XZ];
  Sp[30] += S[np][M_YY];
  Sp[31] += S[np][M_YZ];
  Sp[35] += S[np][M_ZZ];

  Sp[39] += S[np][M_XX];
  Sp[40] += S[np][M_XY];
  Sp[41] += S[np][M_XZ];
  Sp[42] += S[np][M_YY];
  Sp[43] += S[np][M_YZ];
  Sp[44] += S[np][M_ZZ];
}


void ElementSOLID81::make_Omega (uint16 np, Mat<6,9> &B) {
  double *Omega = B.ptr();
  for (uint16 i=0; i < 3; i++) {
    Omega[0*9+(i*3+0)] = O[np][0+i*3]; // exx
    Omega[1*9+(i*3+0)] = O[np][1+i*3]; // exy
    Omega[1*9+(i*3+1)] = O[np][0+i*3]; // exy
    Omega[2*9+(i*3+0)] = O[np][2+i*3]; // exz
    Omega[2*9+(i*3+2)] = O[np][0+i*3]; // exz
    Omega[3*9+(i*3+1)] = O[np][1+i*3]; // eyy
    Omega[4*9+(i*3+1)] = O[np][2+i*3]; // eyz
    Omega[4*9+(i*3+2)] = O[np][1+i*3]; // eyz
    Omega[5*9+(i*3+2)] = O[np][2+i*3]; // ezz
  }
}

bool ElementSOLID81::getScalar(double* scalar, scalarQuery query, uint16 gp, const double scale) {
  //see queries in query.h
  //gp - needed gauss point 
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

  // here we obtain results for particular Gaussian point gp
  assert (gp < nOfIntPoints());

  Vec<3> tmp;
  Mat_Hyper_Isotrop_General* mat;
  double J;
  switch (query) {
    case scalarQuery::SP:
      *scalar += storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE) * scale;
      return true;

    case scalarQuery::WU:
      getVector(&tmp, vectorQuery::IC, gp, 1.0);
      mat = CHECK_NOTNULL(dynamic_cast<Mat_Hyper_Isotrop_General*>(storage->getMaterial()));
      *scalar += mat->W(tmp[0], tmp[1], tmp[2]) * scale;
      return true;

    case scalarQuery::WP:
      J = solidmech::J_C(C[gp].ptr());
      mat = CHECK_NOTNULL(dynamic_cast<Mat_Hyper_Isotrop_General*>(storage->getMaterial()));
      *scalar += 0.5 * mat->getK() * (J - 1.0) * (J - 1.0) * scale;
      return true;

    case scalarQuery::VOL:
      *scalar += volume() * scale;
      return true;
  }
  return false;
}


bool ElementSOLID81::getVector(math::Vec<3>* vector, vectorQuery query, uint16 gp, const double scale) {
  if (gp == GP_MEAN) { //need to average result over the element
    double dWtSum = volume();
    double dWt;
    for (uint16 np = 0; np < nOfIntPoints(); np ++) {
      dWt = intWeight(np);
      bool ret = getVector(vector, query, np, dWt/dWtSum*scale );
      if(ret == false) return false;
    }
    return true;
  }

  // here we obtain results for particular Gaussian point gp
  assert (gp < nOfIntPoints());

  double IC[3];
  switch (query) {
    case vectorQuery::IC:
      solidmech::IC_C(C[gp].ptr(), IC);
      (*vector)[0] += IC[0] * scale;
      (*vector)[1] += IC[1] * scale;
      (*vector)[2] += IC[2] * scale;
      return true;
  }
  return false;
}



//return a tensor in a global coordinate system
bool  ElementSOLID81::getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale) {
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

  assert (gp < nOfIntPoints());

  Mat<3,3> matF;
  MatSym<3> matS;
  double J;
  double cInv[6];
  double pe;

  switch (query) {
    case tensorQuery::COUCHY:

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
      pe = storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE);
      // TODO: it seems that S[gp] contains deviatoric + pressure already..
      // we dont need to sum it again
      for (uint16 i = 0; i < 6; i++) {
        matS.data[i] = S[gp][i] + pe * J * cInv[i];
      }
      matBTDBprod (matF, matS, 1.0/J*scale, *tensor); //Symmetric Couchy tensor
      return true;

    case tensorQuery::PK2:
      // hydrostatic part of S: Sp = p * J * C^(-1)
      J = solidmech::J_C(C[gp].ptr());
      solidmech::invC_C (C[gp].ptr(), J, cInv); 
      pe = storage->getElementDofSolution(getElNum(), Dof::HYDRO_PRESSURE);
      for (uint16 i = 0; i < 6; i++) {
        tensor->data[i] += (S[gp][i] + pe * J * cInv[i]) * scale;
      }
      return true;

    case tensorQuery::C:
      tensor->data[0] += C[gp][M_XX]*scale;
      tensor->data[1] += C[gp][M_XY]*scale;
      tensor->data[2] += C[gp][M_XZ]*scale;
      tensor->data[3] += C[gp][M_YY]*scale;
      tensor->data[4] += C[gp][M_YZ]*scale;
      tensor->data[5] += C[gp][M_ZZ]*scale;
      return true;

    case tensorQuery::E:
      tensor->data[0] += (C[gp][M_XX]-1.0)*0.5*scale;
      tensor->data[1] += C[gp][M_XY]*0.5*scale;
      tensor->data[2] += C[gp][M_XZ]*0.5*scale;
      tensor->data[3] += (C[gp][M_YY]-1.0)*0.5*scale;
      tensor->data[4] += C[gp][M_YZ]*0.5*scale;
      tensor->data[5] += (C[gp][M_ZZ]-1.0)*0.5*scale;
      return true;
  }
  return false;
}

} // namespace nla3d

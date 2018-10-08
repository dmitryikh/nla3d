#include "elements/INTER3.h"

namespace nla3d {

ElementINTER3::ElementINTER3 () {
  type = ElementType::INTER3;
}

void ElementINTER3::pre () {
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY, Dof::UZ});
  }
}

void ElementINTER3::buildK() {
  Eigen::MatrixXd Ke(18,18);
  Ke.setZero();

  Eigen::MatrixXd D(3,3);
  make_D(D);

  makeJacob();

  // build Ke in local sys
  double dWt; //Gaussian quadrature weight
  for (uint16 np=0; np < nOfIntPoints(); np++) {
     for (uint16 npj=0; npj < nOfIntPoints(); npj++) {
      dWt = radoIntWeight(np,npj);
      Eigen::MatrixXd matB = make_B(np,npj);
      Ke.triangularView<Eigen::Upper>() += dWt*matB.transpose()*D*matB;
    }
  }// loop over integration points

  Eigen::MatrixXd T_Ke(18,18);
  Eigen::MatrixXd T_trans = make_T().transpose();
  Eigen::MatrixXd z = Eigen::MatrixXd::Zero(3,3);


  T_Ke << T_trans, z, z, z, z, z,
          z, T_trans, z, z, z, z,
          z, z, T_trans, z, z, z,
          z, z, z, T_trans, z, z,
          z, z, z, z, T_trans, z,
          z, z, z, z, z, T_trans;

  // build Ke in global sys
  Eigen::MatrixXd Ke_glob(18,18);
  Ke_glob.setZero();
  Ke_glob =  T_Ke.transpose()*Ke*T_Ke;     
  assembleK(Ke_glob, {Dof::UX, Dof::UY, Dof::UZ});
}

void ElementINTER3::update () {
  //Перемещения в узлах верхнего и нижнего треугольника
  Eigen::VectorXd U1(9);
  Eigen::VectorXd U2(9);
  U1.setZero();
  U2.setZero();
  for (uint16 i = 0; i < 3; i++) {
    U1(i*3 + 0) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U1(i*3 + 1) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
    U1(i*3 + 2) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UZ);
  }
  for (uint16 i = 0; i < 3; i++) {
    U2(i*3 + 0) = storage->getNodeDofSolution(getNodeNumber(i+3), Dof::UX);
    U2(i*3 + 1) = storage->getNodeDofSolution(getNodeNumber(i+3), Dof::UY);
    U2(i*3 + 2) = storage->getNodeDofSolution(getNodeNumber(i+3), Dof::UZ);
  }

  Eigen::MatrixXd T = make_T();
  Eigen::MatrixXd T_inv = T.inverse();
  Eigen::MatrixXd T_U(9,9);
  Eigen::MatrixXd z = Eigen::MatrixXd::Zero(3,3);

  T_U << T_inv, z, z,
         z, T_inv, z,
         z, z, T_inv;
  U1 = T_U*U1;
  U2 = T_U*U2;

  //Интегрированные по площади перемещения верхнего и нижнего треугольников
  Eigen::VectorXd U1S(3);
  Eigen::VectorXd U2S(3);
  U1S.setZero();
  U2S.setZero();

  for (uint16 np=0; np < nOfIntPoints(); np++) {
    for (uint16 npj=0; npj < nOfIntPoints(); npj++) {
      double dWt = radoIntWeight(np,npj);
      Eigen::MatrixXd matB = make_subB(np,npj);
      U1S += dWt*matB*U1;
      U2S += dWt*matB*U2;
    }
  }

  Eigen::VectorXd strainsE(3);
  strainsE = U2S - U1S;

  Eigen::MatrixXd D(3,3);
  make_D(D);
  Eigen::VectorXd stressE(3);
  stressE = D*strainsE;

  // переход в глобальную ск
  strainsE = T*strainsE;
  stressE = T*stressE;

  for (int i(0); i < 3; i++){
    stress[i] = stressE(i);
    strains[i] = strainsE(i);
  }
}

void ElementINTER3::makeJacob(){
  //Находим координаты узлов в локальной декартовой системе координат треугольника (t1, t2)
  //Начало координт - первый узел треугольника
  math::Vec<2> locX1(0.,0.);
  //Вектор s1 направлен по одной из сторон
  math::Vec<3> t1 = storage->getNode(getNodeNumber(1)).pos-storage->getNode(getNodeNumber(0)).pos;
  math::Vec<3> t2 = storage->getNode(getNodeNumber(2)).pos-storage->getNode(getNodeNumber(0)).pos;

  math::Vec<2> locX2(t1.length(),0.);
  //Для координат третьего узла нужна ортогонолизация. Ищем угол треугольника при начале координат
  double mult = t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2];
  double angle = acos(mult/t1.length()/t2.length());
  math::Vec<2> locX3(t2.length()*cos(angle),t2.length()*sin(angle));

  math::Mat<2,2> J(locX1[0]-locX3[0],locX1[1]-locX3[1],
                   locX2[0]-locX3[0],locX2[1]-locX3[1]);

  det = J.det(); //Якобиан перехода между L координатами и локальными декартовыми
  if (det < 0)
    LOG(FATAL) << "Negative Jacobian value " << det;
}

void ElementINTER3::make_D(Eigen::MatrixXd& D){
  D << ks, 0., 0.,
       0., ks, 0., 
       0., 0., kn;
}

Eigen::MatrixXd ElementINTER3::make_subB(uint16 np, uint16 npj){

  Eigen::MatrixXd N1(3,3);
  Eigen::MatrixXd N2(3,3);
  Eigen::MatrixXd N3(3,3);

  N1 << radoIntL1(np), 0., 0., 
        0., radoIntL1(np), 0., 
        0., 0., radoIntL1(np);
  N2 << radoIntL2(np,npj), 0., 0., 
        0., radoIntL2(np,npj), 0., 
        0., 0., radoIntL2(np,npj);
  N3 << radoIntL3(np,npj), 0., 0., 
        0., radoIntL3(np,npj), 0., 
        0., 0., radoIntL3(np,npj);    

  Eigen::MatrixXd B(3,9);
  B << N1, N2, N3;

 return B;

}

Eigen::MatrixXd ElementINTER3::make_B(uint16 np, uint16 npj){ 
  // in local sys
  Eigen::MatrixXd N1(3,3);
  Eigen::MatrixXd N2(3,3);
  Eigen::MatrixXd N3(3,3);

  N1 << radoIntL1(np), 0., 0., 
        0., radoIntL1(np), 0., 
        0., 0., radoIntL1(np);
  N2 << radoIntL2(np,npj), 0., 0., 
        0., radoIntL2(np,npj), 0., 
        0., 0., radoIntL2(np,npj);
  N3 << radoIntL3(np,npj), 0., 0., 
        0., radoIntL3(np,npj), 0., 
        0., 0., radoIntL3(np,npj);     

  Eigen::MatrixXd B(3,18);
  B << N1, N2, N3, -N1, -N2, -N3;

  return B;
}

Eigen::MatrixXd ElementINTER3::make_T(){
  //Востанавливаем локальный базис s1,s2,n
  //s1 совпадает с одной из сторон
  math::Vec<3> s1 = storage->getNode(getNodeNumber(1)).pos - storage->getNode(getNodeNumber(0)).pos;
  math::Vec<3> t2 = storage->getNode(getNodeNumber(2)).pos - storage->getNode(getNodeNumber(0)).pos;
  //Востанавливаем нормаль как векторное произведение двух вектров в плоскости треугольника
  math::Vec<3> n(s1[1]*t2[2]-s1[2]*t2[1],s1[2]*t2[0]-s1[0]*t2[2],s1[0]*t2[1]-s1[1]*t2[0]);
  //Востанавливаем s2 как векторное произведение n x s1
  math::Vec<3> s2(n[1]*s1[2]-n[2]*s1[1],n[2]*s1[0]-n[0]*s1[2],n[0]*s1[1]-n[1]*s1[0]);

  n = n*(1./n.length());
  s1 = s1*(1./s1.length());
  s2 = s2*(1./s2.length());
  
  //Матрица поворота от локальной к глобальной ск
  Eigen::MatrixXd T(3,3); 
  T << s1[0],s2[0],n[0],
       s1[1],s2[1],n[1],
       s1[2],s2[2],n[2];

  normal = n;

  return T;
}

bool  ElementINTER3::getVector(math::Vec<3>* vector, vectorQuery query, uint16 gp, const double scale) {
  return false;
}

bool ElementINTER3::getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale){
  if (query == tensorQuery::C){
      tensor->comp(0,0) += strains[0];
      tensor->comp(1,1) += strains[1];
      tensor->comp(2,2) += strains[2];
      tensor->comp(0,1) += 0.;
      tensor->comp(1,2) += 0.;
      tensor->comp(0,2) += 0.;
      return true;
  }
  if (query == tensorQuery::E){
    tensor->comp(0,0) += stress[0];
    tensor->comp(1,1) += stress[1];
    tensor->comp(2,2) += stress[2];
    tensor->comp(0,1) += 0.;
    tensor->comp(1,2) += 0.;
    tensor->comp(0,2) += 0.;
    return true;
  }
}


} //namespace nla3d

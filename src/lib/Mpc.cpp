// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "Mpc.h"
#include "FEStorage.h"
#include "Node.h"

namespace nla3d {

using namespace math;


void Mpc::print (std::ostream& out) {
  auto token = eq.begin();
  bool first = true;
  while (token != eq.end()) {
    //Cij_add(eq_num, token->node, token->node_dof, token->coef);
    if (!first) {
      out << " + ";
    } else {
      first = false;
    }
    out << toStr(token->coef) << " * " << toStr(token->node) << ":" << Dof::dofTypeLabels[token->node_dof];
    token++;
  }
  out << " = " << toStr(b);
}

void MpcCollection::registerMpcsInStorage () {
  for (size_t i = 0; i < collection.size(); i++) {
    storage->addMpc(collection[i]);
  }
}


void MpcCollection::printEquations (std::ostream& out) {
  for (size_t i = 0; i < collection.size(); i++) {
    collection[i]->print(out);
    out << std::endl;
  }
}

RigidBodyMpc::RigidBodyMpc () {

}

RigidBodyMpc::~RigidBodyMpc () {
  collection.clear();
}

void RigidBodyMpc::pre () {
  // create Mpc for all nodes and all dofs..
  collection.assign(slaveNodes.size() * dofs.size(), NULL);
  for (size_t i = 0; i < slaveNodes.size(); i++) {
    for (size_t j = 0; j < dofs.size(); j++) {
      Mpc* mpc = new Mpc;
      mpc->b = 0.0;
      mpc->eq.push_back(MpcTerm(slaveNodes[i], dofs[j], 1.0));
      mpc->eq.push_back(MpcTerm(masterNode, dofs[j], -1.0));
      mpc->eq.push_back(MpcTerm(masterNode, Dof::ROTX, 0.0));
      mpc->eq.push_back(MpcTerm(masterNode, Dof::ROTY, 0.0));
      mpc->eq.push_back(MpcTerm(masterNode, Dof::ROTZ, 0.0));
      collection[i*dofs.size() + j] = mpc;
    }
  }
  storage->addNodeDof(masterNode, {Dof::UX, Dof::UY, Dof::UZ, Dof::ROTX, Dof::ROTY, Dof::ROTZ});
}

void RigidBodyMpc::update () {
  // rigid body motion formula is the next:
  // {u} = {w} +([C]-[I])*({p}-{q})
  // where  {u} - displacement of a slave node
  //        {w} - displacement of the master node
  //        [C] - rotation matrix (depends on rotation vector {theta})
  //        {p} - position of the master node
  //        {q} - position of a slave node
  // where [C] in componentwise form:
  // C_{ij} = \delta_{ij} cos\theta + \frac{sin\theta}{\theta}e_{ikj}\theta_k + \left(\frac{1-cos\theta}{\theta^2}\theta_i \theta_j\right) 
  Vec<3> theta0;
  Vec<3> w0;
  theta0[0] = storage->getNodeDofSolution(masterNode, Dof::ROTX);
  theta0[1] = storage->getNodeDofSolution(masterNode, Dof::ROTY);
  theta0[2] = storage->getNodeDofSolution(masterNode, Dof::ROTZ);
  w0[0] = storage->getNodeDofSolution(masterNode, Dof::UX);
  w0[1] = storage->getNodeDofSolution(masterNode, Dof::UY);
  w0[2] = storage->getNodeDofSolution(masterNode, Dof::UZ);
  double thetaNorm = theta0.length();
  Vec<3> masterPos;
  Vec<3> slavePos;
  Mat<3,3> C;
  double c1, c2, c3, c4;
  if (thetaNorm < 1.0e-3) {
    // for very small angles we need to define some undefined values
    // TODO: for now eps = 1.0e-3 was choosed randomly, need to check it
    c1 = 1.0;
    c2 = 0.5;
    c3 = -1.0/3.0;
    c4 = -1.0/12.0;
  } else {
  // define constant exactly
    c1 = sin(thetaNorm) / thetaNorm;
    c2 = (1.0 - cos(thetaNorm)) / (thetaNorm * thetaNorm);
    c3 = (thetaNorm * cos(thetaNorm) - sin(thetaNorm)) / (thetaNorm * thetaNorm * thetaNorm);
    c4 = (thetaNorm * sin(thetaNorm) - 2.0 + 2.0 * cos(thetaNorm)) / (thetaNorm * thetaNorm * thetaNorm * thetaNorm);
  }
  for (uint16 i = 0; i < 3; i++) {
    for (uint16 j = 0; j < 3; j++) {
      double theta0Mat_ij = solidmech::LeviCivita[i][0][j] * theta0[0] +
                      solidmech::LeviCivita[i][1][j] * theta0[1] +
                      solidmech::LeviCivita[i][2][j] * theta0[2];
      C[i][j] = cos(thetaNorm) * solidmech::I[i][j] + c1 * theta0Mat_ij + 
              c2 * theta0[i] * theta0[j];
    }
  }
  storage->getNodePosition(masterNode, masterPos.ptr(), false);
  for (size_t n = 0; n < slaveNodes.size(); n++) {
    storage->getNodePosition(slaveNodes[n], slavePos.ptr(), false);
    Vec<3> pqVec = slavePos - masterPos;
    // drivatives for C matrix for i row (with respect to j and l)
    Mat<3,3> cDeriv;
    for (size_t d = 0; d < dofs.size(); d++) {
      uint16 i;
      switch (dofs[d]) {
        case Dof::UX:
          i = 0;
          break;
        case Dof::UY:
          i = 1;
          break;
        case Dof::UZ:
          i = 2;
          break;
        default:
          LOG(FATAL) << "which degree of freedom?";
        }
        for (uint16 j = 0; j < 3; j++) {
          for (uint16 l = 0; l < 3; l++) {
            cDeriv[j][l] = solidmech::I[i][j]*(-c1)*theta0[l] + c3*theta0[l]*(
                  solidmech::LeviCivita[i][0][j]*theta0[0] +
                  solidmech::LeviCivita[i][1][j]*theta0[1] +
                  solidmech::LeviCivita[i][2][j]*theta0[2]) +
                c1*solidmech::LeviCivita[i][l][j] +
                c4*theta0[l]*theta0[i]*theta0[j] +
                c2*(solidmech::I[i][l]*theta0[j] + theta0[i]*solidmech::I[j][l]);
          }
        }
        collection[n*static_cast<uint16> (dofs.size()) + d]->b = w0[i] + (C[i][0] - solidmech::I[i][0])*pqVec[0] +
              (C[i][1] - solidmech::I[i][1])*pqVec[1] + (C[i][2] - solidmech::I[i][2])*pqVec[2] -
              storage->getNodeDofSolution(slaveNodes[n], dofs[d]);
        auto token = collection[n*dofs.size() + d]->eq.begin();
        token++;
        token++;
        // masterNode::ROTX coef
        token->coef = - (cDeriv[0][0]*pqVec[0] + cDeriv[1][0]*pqVec[1] + cDeriv[2][0]*pqVec[2]);
        token++;
        // masterNode::ROTY coef
        token->coef = - (cDeriv[0][1]*pqVec[0] + cDeriv[1][1]*pqVec[1] + cDeriv[2][1]*pqVec[2]);
        token++;
        // masterNode::ROTZ coef
        token->coef = - (cDeriv[0][2]*pqVec[0] + cDeriv[1][2]*pqVec[1] + cDeriv[2][2]*pqVec[2]);
    }
  }
}


} // namespace nla3d

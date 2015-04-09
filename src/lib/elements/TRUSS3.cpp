#include "elements/TRUSS3.h"

namespace nla3d {

ElementTRUSS3::ElementTRUSS3 () {
  Element::number_of_nodes = 2;
  number_of_dimensions = 3;
  nodes = new uint32[Element::number_of_nodes];

  Node::registerDofType(Dof::UX);
  Node::registerDofType(Dof::UY);
  Node::registerDofType(Dof::UZ);
}

void ElementTRUSS3::pre () {
  for (uint16 i = 0; i < Element::n_nodes(); i++) {
    storage->registerNodeDof(getNodeNumber(i), Dof::UX);
    storage->registerNodeDof(getNodeNumber(i), Dof::UY);
    storage->registerNodeDof(getNodeNumber(i), Dof::UZ);
  }
}

// here stiffness matrix is building with Eignen's library matrix routines
void ElementTRUSS3::build () {
  // Ke will store element stiffness matrix
  Eigen::MatrixXd Ke(6,6);
  // T for transformation matrix (to map governing equations from an element local coordinate system
  // into global coordinate system
  Eigen::MatrixXd T(2,6);
  // K stores element stiffness matrix in a local coordinate system
  Eigen::MatrixXd K(2,2);

  Ke.setZero();
  T.setZero();
  K.setZero();

  // here we use pointer to FEStorage class where all FE information is stored (from node and
  // element tables to solution results). getNode(..).pos used to obtain an class instance for
  // particular node an get its position which is Vec<3> data type. 
  math::Vec<3> deltaPos = storage->getNode(getNodeNumber(1)).pos - storage->getNode(getNodeNumber(0)).pos;
  // as long as we need 1.0/length instead of length itself it's useful to store this value in the
  // variable inv_length. Such kind of code optimization is quite useful in massively computational
  // expensive procedures.
  double inv_length = 1.0 / deltaPos.length();

  // filling transformation matrix T (see reference theretical material to get it clear. The link to
  // theoretical material could be found in TRUSS3.h file). Here Eigen library matrix and its
  // operator(). For more information read Eigen manual.
  T(0,0) = deltaPos[0] * inv_length;
  T(0,1) = deltaPos[1] * inv_length;
  T(0,2) = deltaPos[2] * inv_length;
  T(1,3) = T(0,0);
  T(1,4) = T(0,1);
  T(1,5) = T(0,2);

  // This is just another way in Eigen to initialize {{1.0, -1.0}, {-1.0, 1.0}} matrix. 
  K << 1.0, -1.0,
      -1.0,  1.0;
  K *= A*E*inv_length;
  // Transform local stiffness matrix K into the global coordinate sytem (matrix Ke). As far as Ke
  // is symmetric it's ok to calculate just upper triangular part of it. For this reason here
  // Ke.triangularView<Eigen::Upper> is used.
  Ke.triangularView<Eigen::Upper>() = T.transpose() * K * T;

  // start assemble procedure. Here we should provide element stiffness matrix and an order of 
  // nodal DoFs in the matrix.
  assemble(Ke, {Dof::UX, Dof::UY, Dof::UZ});
}

// after solution it's handy to calculate stresses, strains and other stuff in elements. In this
// case a truss stress will be restored
void ElementTRUSS3::update () {
  // T for transformation matrix
  Eigen::MatrixXd T(2,6);
  // B for strain matrix
  Eigen::MatrixXd B(1,2);

  T.setZero();
  B.setZero();

  math::Vec<3> deltaPos = storage->getNode(getNodeNumber(1)).pos - storage->getNode(getNodeNumber(0)).pos;
  double inv_length = 1.0 / deltaPos.length();

  T(0,0) = deltaPos[0] * inv_length;
  T(0,1) = deltaPos[1] * inv_length;
  T(0,2) = deltaPos[2] * inv_length;
  T(1,3) = T(0,0);
  T(1,4) = T(0,1);
  T(1,5) = T(0,2);

  B << -inv_length, inv_length;

  // read solution results from FEStorage. Here we need to fill nodal DoFs values into U vector.
  Eigen::VectorXd U(6);
  for (uint16 i = 0; i < Element::n_nodes(); i++) {
    U(i*3 + 0) = storage->getDofSolution(getNodeNumber(i), Dof::UX);
    U(i*3 + 1) = storage->getDofSolution(getNodeNumber(i), Dof::UY);
    U(i*3 + 2) = storage->getDofSolution(getNodeNumber(i), Dof::UZ);
  }

  // The formula to calculate truss stresses (see theory reference):
  // S = E  *  B  *  T  *  D
  //1x1=1x1 * 1x2 * 2x6 * 6x1
  S = E * (B * T * U)(0,0);
}

} //namespace nla3d

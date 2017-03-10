// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/element.h"

namespace nla3d {


Element::Element () {

}


Element::~Element() {
  if (nodes) {
    delete[] nodes;
    nodes = nullptr;
  }
}


void Element::print (std::ostream& out) {
  out << "E " << getElNum() << ":";
  for (uint16 i = 0; i < getNNodes(); i++) {
    out << "\t" << getNodeNumber(i);
  }
}

Element& Element::operator= (const Element& from)
{
  assert(nodes && from.nodes);
  memcpy(nodes, from.nodes, sizeof(uint32)*getNNodes());
  return *this;
}


void Element::assembleK(Eigen::Ref<Eigen::MatrixXd> Ke,
                       std::initializer_list<Dof::dofType> _nodeDofs) {
  assert (nodes != NULL);
  assert (Ke.rows() == Ke.cols());
  std::vector<Dof::dofType> nodeDof(_nodeDofs);
  uint16 dim = static_cast<uint16> (_nodeDofs.size());
  assert (getNNodes() * dim == Ke.rows());

  for (uint16 i=0; i < getNNodes(); i++) {
    for (uint16 di=0; di < dim; di++) {
      for (uint16 j=i; j < getNNodes(); j++) {
        for (uint16 dj=0; dj < dim; dj++) {
          if ((i==j) && (dj<di)) {
            continue;
          } else {
            storage->addValueK(nodes[i], nodeDof[di], nodes[j], nodeDof[dj], 
                Ke.selfadjointView<Eigen::Upper>()(i*dim+di, j*dim +dj));
          }
        }
      }
    }
  }
}


void Element::buildC() {
  LOG(FATAL) << "buildC is not implemented";
}


void Element::buildM() {
  LOG(FATAL) << "buildM is not implemented";
}

bool Element::getScalar(double* scalar, scalarQuery code, uint16 gp, const double scale) {
  // TODO: check that LOG_N_TIMES macro work correctly inside of virtual functions
  LOG_N_TIMES(10, WARNING) << "getScalar function is not implemented";
  return false;
}

bool Element::getVector(math::Vec<3>* vector, vectorQuery code, uint16 gp, const double scale) {
  LOG_N_TIMES(10, WARNING) << "getVector function is not implemented";
  return false;
}

bool  Element::getTensor(math::MatSym<3>* tensor, tensorQuery code, uint16 gp, const double scale) {
  LOG_N_TIMES(10, WARNING) << "getTensor function is not implemented";
  return false;
}

} // namespace nla3d

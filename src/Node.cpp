// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "Node.h"

namespace nla3d {

//TODO: could it be used special cosntructor for vector to define right here?
std::vector<uint16> Node::dofNumberList;
uint16 Node::numberOfDofs = 0;

Node::Node() {
  // initialise static variables (just once)
  if (dofNumberList.size() == 0) {
    dofNumberList.assign(Dof::numberOfDofTypes, Dof::UNDEFINED);
    numberOfDofs = 0;
  }
}

void Node::registerDofType(Dof::dofType type) {
  bool isFound = false;
  if (dofNumberList[type] == Dof::UNDEFINED) {
    dofNumberList[type] = numberOfDofs;
    numberOfDofs++;
  }
}

Dof::dofType Node::getDofType (uint16 dofIndex) {
  assert(dofIndex < getNumberOfDofs());
  for (size_t i = 0; i < Dof::numberOfDofTypes; i++) {
    if (getDofIndex((Dof::dofType)i) == dofIndex) {
      return (Dof::dofType) i;
    }
  }
  error ("What I'm doing here?");
  return Dof::UNDEFINED;
}

}

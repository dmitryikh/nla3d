// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#include "Node.h"

namespace nla3d {


void Node::display (uint32 nn) {
  LOG(INFO) << "N " << nn << " : " << pos;
}


std::string Node::toString() {
  std::string str;
  char buff[100];
  sprintf_s(buff,100,"%f %f %f",pos[0], pos[1], pos[2]);
  str+=buff;
  return str; //TODO: do it easy
}

}

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"

namespace nla3d {
class Element;

class ElementFactory {
  public:
  enum elTypes {
    NOT_DEFINED,
    PLANE41,
    SOLID81,
    TRUSS3,
    LAST
  };

  static const char* const elTypeLabels[];

  static elTypes elName2elType (std::string elName); 
  static void createElements (elTypes elId, const uint32 n, std::vector<Element*>& ptr); 
};

} // namespace nla3d

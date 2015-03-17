// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"

namespace nla3d {

class ElementFactory {
  public:
  enum elTypes {
    NOT_DEFINED,
    PLANE41,
    SOLID81,
    LAST
  };

  static const char* const elTypeLabels[];

  static elTypes elName2elType (std::string elName); 
  static void createElements (elTypes elId, const uint32 n, std::vector<Element*>& ptr); 
};

} // namespace nla3d

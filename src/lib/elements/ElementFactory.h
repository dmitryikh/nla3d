// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"

namespace nla3d {
class Element;

class ElementFactory {
  public:

    static ElementType elName2elType (std::string elName); 
    static void createElements (ElementType elId, const uint32 n, std::vector<Element*>& ptr); 
};

} // namespace nla3d

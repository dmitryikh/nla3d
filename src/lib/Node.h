// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#pragma once
#include "sys.h"
#include "math/Vec.h"
#include "Dof.h"

namespace nla3d {

//class Node represents spatial 3D node
class Node {
  public:
    //in-out operation: (rudiment actually..)
    void display (uint32 nn);
    std::string toString();

    math::Vec<3> pos;

    friend class Element;
};

} // namespace nla3d

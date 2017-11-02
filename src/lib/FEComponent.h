// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#pragma once
#include "sys.h"

namespace nla3d {

class FEComponent {
  public:
    enum typeOfComponent {
      NOT_DEFINED,
      ELEMENTS,
      NODES,
      LAST
    };

    static const char* const labelsOfComponent[];

    FEComponent ();
    ~FEComponent ();

    std::vector<uint32> list;
    typeOfComponent type;
    std::string name;

    static typeOfComponent typeFromString(const std::string& typeName);
    void print ();
};

inline MAKE_LOGGABLE(FEComponent, obj, os) {
  os << obj.name << ": " << obj.list.size() << " " << obj.labelsOfComponent[obj.type];
  return os;
}

} // namespace nla3d

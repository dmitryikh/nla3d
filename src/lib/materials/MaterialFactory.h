// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#pragma once
#include "materials/materials_hyperelastic.h"

namespace nla3d {

class MaterialFactory {
  public:
  enum matId {
    NOT_DEFINED,
    NEO_HOOKEAN_COMP,
    BIDERMAN_COMP,
    MOONEYRIVLIN_COMP,
    LAST
  };

  static const char* const matModelLabels[];

  static matId matName2matId (std::string matName);
  static Material* createMaterial (std::string matName);
};

} // namespace nla3d

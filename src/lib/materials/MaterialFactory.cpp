// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#include "materials/MaterialFactory.h"

namespace nla3d {

const char* const MaterialFactory::matModelLabels[]={"UNDEFINED",
  "Neo-Hookean",
  "Biderman",
  "Mooney-Rivlin"
};

MaterialFactory::matId MaterialFactory::matName2matId (std::string matName) {
  for (uint16 i = 0; i < MaterialFactory::LAST; i++) {
    if (matName.compare(MaterialFactory::matModelLabels[i]) == 0) {
      return (MaterialFactory::matId) i;
    }
  }
  return MaterialFactory::NOT_DEFINED;
}

Material* MaterialFactory::createMaterial (std::string matName) {
  uint16 matId = matName2matId(matName);
  Material* mat;
  if (matId == MaterialFactory::NOT_DEFINED) {
    LOG(ERROR) << "Can't find material " << matName;
    return nullptr;
  }
  switch (matId) {
    case MaterialFactory::NEO_HOOKEAN_COMP:
      mat = new Mat_Comp_Neo_Hookean();
      break;
    case MaterialFactory::BIDERMAN_COMP:
      mat = new Mat_Comp_Biderman();
      break;
    case MaterialFactory::MOONEYRIVLIN_COMP:
      mat = new Mat_Comp_MooneyRivlin();
      break;
    default:
      LOG(ERROR) << "Don't have a material with id  = " << matId;
      return nullptr;
  }
  mat->code = matId;
  return mat;
}

} // namespace nla3d

#include "materials/material_factory.h"

const char* const MaterialFactory::matModelLabels[]={"UNDEFINED",
  "Neo-Hookean",
  "Biderman",
  "Mooney-Rivlin"
};

MaterialFactory::matId MaterialFactory::matName2matId (string matName) {
  for (uint16 i = 0; i < MaterialFactory::LAST; i++) {
    if (matName.compare(MaterialFactory::matModelLabels[i]) == 0) {
      return (MaterialFactory::matId) i;
    }
  }
  return MaterialFactory::NOT_DEFINED;
}

Material* MaterialFactory::createMaterial (string matName) {
  uint16 matId = matName2matId(matName);
  Material* mat;
  if (matId == MaterialFactory::NOT_DEFINED)
    error("createMaterial: can't find material %s", matName.c_str());
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
      error("createMaterial: don't have a material with id %d", matId);
  }
  mat->code = matId; //TODO: its could be a bad thing to pass matId here
  return mat;
}


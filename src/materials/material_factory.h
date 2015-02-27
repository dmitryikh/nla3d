#pragma once
#include "materials/materials_hyperelastic.h"

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

  static matId matName2matId (string matName); 
  static Material* createMaterial (string matName); 
};

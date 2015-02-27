#pragma once
#include "elements/element.h"

class ElementFactory {
  public:
  enum elTypes {
    NOT_DEFINED,
    PLANE41,
    SOLID81,
    LAST
  };

  static const char* const elTypeLabels[];

  static elTypes elName2elType (string elName); 
  static void createElements (elTypes elId, const uint32 n, vector<Element*>& ptr); 
};

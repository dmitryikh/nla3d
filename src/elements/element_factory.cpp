#include "elements/element_factory.h"
#include "elements/PLANE41.h"
#include "elements/SOLID81.h"

const char* const ElementFactory::elTypeLabels[]={"UNDEFINED",
  "PLANE41",
  "SOLID81",
};

ElementFactory::elTypes ElementFactory::elName2elType (string elName) {
  for (uint16 i = 1; i < (uint16) ElementFactory::LAST; i++) {
    if (elName.compare(ElementFactory::elTypeLabels[i]) == 0) {
      return (ElementFactory::elTypes) i;
    }
  }
  return ElementFactory::NOT_DEFINED;
}
void ElementFactory::createElements (elTypes elId, const uint32 n, vector<Element*>& ptr) {
  if (elId == ElementFactory::NOT_DEFINED)
    error("ElementFactory::createElements: element type %s is undefined");
  switch (elId) {
    case ElementFactory::PLANE41:
      for (uint32 i = 0; i < n; i++) {
        ptr.push_back(new ElementPLANE41());
      }
      break;
    case ElementFactory::SOLID81:
      for (uint32 i = 0; i < n; i++) {
        ptr.push_back(new ElementSOLID81());
      }
      break;
    default:
      error("ElementFactory::createElements: don't have an element with id %d", elId);
  }
}



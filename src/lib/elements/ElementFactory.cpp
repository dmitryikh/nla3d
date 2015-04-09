// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/ElementFactory.h"
#include "elements/PLANE41.h"
#include "elements/SOLID81.h"
#include "elements/TRUSS3.h"

namespace nla3d {

const char* const ElementFactory::elTypeLabels[]={"UNDEFINED",
  "PLANE41",
  "SOLID81",
  "TRUSS3"
};

ElementFactory::elTypes ElementFactory::elName2elType (std::string elName) {
  for (uint16 i = 1; i < (uint16) ElementFactory::LAST; i++) {
    if (elName.compare(ElementFactory::elTypeLabels[i]) == 0) {
      return (ElementFactory::elTypes) i;
    }
  }
  return ElementFactory::NOT_DEFINED;
}
void ElementFactory::createElements (elTypes elId, const uint32 n, std::vector<Element*>& ptr) {
  if (elId == ElementFactory::NOT_DEFINED) {
    LOG(WARNING) << "Element type is undefined";
  }
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
    case ElementFactory::TRUSS3:
      for (uint32 i = 0; i < n; i++) {
        ptr.push_back(new ElementTRUSS3());
      }
      break;
    default:
      LOG(ERROR) << "Don't have an element with id " << elId;
  }
}

} // namespace nla3d

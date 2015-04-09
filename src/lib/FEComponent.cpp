// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "FEComponent.h"

namespace nla3d {

const char* const FEComponent::labelsOfComponent[]={"NOT DEFINED","ELEM", "NODE", "LAST"};

FEComponent::FEComponent () {
  type = NOT_DEFINED;
  name = "";
}
FEComponent::~FEComponent () {
  list.clear();
} 

FEComponent::typeOfComponent FEComponent::typeFromString(const std::string& typeName) {
  for (size_t i = 1; i < LAST; i++) {
    if (typeName.compare(labelsOfComponent[i]) == 0) {
      return static_cast<typeOfComponent> (i);
    }
  }
  LOG(ERROR) << "Can't find FEComponent::type by name " << typeName;
  return NOT_DEFINED;
}

} // namespace nla3d

#include "Dof.h"

namespace nla3d {

const uint16 Dof::numberOfDofTypes = 7;
const char* const Dof::dofTypeLabels[] = {"UX", "UY", "UZ", "ROTX", "ROTY", 
  "ROTZ", "HYDRO_PRESSURE", "UNDEFINED"};


Dof::dofType Dof::label2dofType (const std::string& label) {
  for (uint16 i = 0; i < numberOfDofTypes; i++) {
    if (label.compare(dofTypeLabels[i]) == 0) {
      return (dofType) i;
    }
  }
  error ("Dof::label2dofType: unknown dof key: %s", label.c_str());
  return Dof::UNDEFINED;
}

} // namespace nla3d

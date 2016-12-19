#include "Dof.h"

namespace nla3d {

const uint16 Dof::numberOfDofTypes = Dof::UNDEFINED;
const char* const Dof::dofTypeLabels[] = {"UX", "UY", "UZ", "ROTX", "ROTY",
  "ROTZ", "HYDRO_PRESSURE", "UNDEFINED"};

// Check that TypeLabels has the same length as numberOfDofTypes
static_assert(Dof::numberOfDofTypes == sizeof(Dof::dofTypeLabels)/sizeof(Dof::dofTypeLabels[0]) - 1);

Dof::dofType Dof::label2dofType (const std::string& label) {
  for (uint16 i = 0; i < numberOfDofTypes; i++) {
    if (label.compare(dofTypeLabels[i]) == 0) {
      return (dofType) i;
    }
  }
  LOG(ERROR) << "Unknown dof key = " << label;
  return Dof::UNDEFINED;
}


void DofCollection::buildDofTable(uint32 _numberOfEntities) {
  //TODO: make it possible to buildDofTable many times
  assert(dofTable.size() == 0);

  numberOfEntities = _numberOfEntities;
  //TODO: not very good idea to initialize all ever possible Dof types at the start..
  dofTable.assign(numberOfEntities * Dof::numberOfDofTypes, Dof());
  numberOfAllocatedDofs = numberOfEntities * Dof::numberOfDofTypes;
  numberOfUsedDofs = 0;
}


// n index starts from 1
void DofCollection::useDof(uint32 n, Dof::dofType dof) {
  assert(n <= numberOfEntities);
  assert(dofTable.size() > 0);
  uint32 ind = (n - 1) * Dof::numberOfDofTypes + dof;
  if (dofTable[ind]._isUsed == false) {
    dofTable[ind]._isUsed = true;
    numberOfUsedDofs++;
  }
}


void DofCollection::clearDofTable() {
  dofTable.clear();
  numberOfUsedDofs = 0;
  numberOfAllocatedDofs = 0;
  numberOfEntities = 0;
}

} // namespace nla3d

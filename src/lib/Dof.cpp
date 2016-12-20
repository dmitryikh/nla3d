#include "Dof.h"

namespace nla3d {

const uint16 Dof::numberOfDofTypes = Dof::UNDEFINED;
const char* const Dof::dofTypeLabels[] = {"UX", "UY", "UZ", "ROTX", "ROTY",
  "ROTZ", "HYDRO_PRESSURE", "UNDEFINED"};

uint32 DofCollection::empty = 0xFFFFFFFF;

// Check that TypeLabels has the same length as numberOfDofTypes
static_assert(Dof::numberOfDofTypes == sizeof(Dof::dofTypeLabels)/sizeof(Dof::dofTypeLabels[0]) - 1,
                "dofTypeLabels and dofType must have the same number of elements");

Dof::dofType Dof::label2dofType (const std::string& label) {
  for (uint16 i = 0; i < numberOfDofTypes; i++) {
    if (label.compare(dofTypeLabels[i]) == 0) {
      return (dofType) i;
    }
  }
  LOG(ERROR) << "Unknown dof key = " << label;
  return Dof::UNDEFINED;
}


void DofCollection::initDofTable(uint32 _numberOfEntities) {
  //TODO: make it possible to buildDofTable many times
  assert(dofPos.size() == 0);
  assert(dofs.size() == 0);

  numberOfEntities = _numberOfEntities;
  // dofPos.assign(numberOfEntities + 1, empty);
  dofPos.assign(numberOfEntities + 1, 0);
  numberOfUsedDofs = 0;
}


// n index starts from 1
//
void DofCollection::addDof(uint32 n, std::initializer_list<Dof::dofType> _dofs) {
  assert(n <= numberOfEntities);
  assert(dofPos.size() > 0);

  std::vector<Dof> newDofs;
  // if already registered - just return
  if (dofPos[n-1] == dofPos[n]) {
    for (auto v : _dofs)
      newDofs.push_back(Dof(v));
  } else {
    for (auto v : _dofs) {
      bool isFound = false;
      for (auto it = dofs.begin() + dofPos[n-1]; it < dofs.begin() + dofPos[n]; it++) {
        if (it->type == v) {
          isFound = true;
          break;
        }
      }
      if (!isFound) {
        newDofs.push_back(Dof(v));
      }
    }
  }

  if (newDofs.size() == 0) return;

  dofs.insert(dofs.begin() + dofPos[n], newDofs.begin(), newDofs.end());
  uint32 i = n;

  // incement dofPos with newDofs.size() for all above
  while (i <= numberOfEntities) dofPos[i++] += newDofs.size();

  numberOfUsedDofs += newDofs.size();
}


void DofCollection::clearDofTable() {
  dofs.clear();
  dofPos.clear();
  numberOfUsedDofs = 0;
  numberOfEntities = 0;
}


// n index starts from 1
bool DofCollection::isDofUsed(uint32 n, Dof::dofType dof) {
  assert(n <= numberOfEntities);
  assert(dofPos.size() > 0);
  if (dofPos[n-1] == empty || dofPos[n] == empty || dofPos[n-1] == dofPos[n]) return false;
  for (auto it = dofs.begin() + dofPos[n-1]; it < dofs.begin() + dofPos[n]; it++) {
    if (it->type == dof) return true;
  }
  return false;
}


Dof* DofCollection::getDof(uint32 n, Dof::dofType dof) {
  assert(n <= numberOfEntities);
  assert(dofPos.size() > 0);

  Dof* pdof = nullptr;

  for (auto it = dofs.begin() + dofPos[n-1]; it < dofs.begin() + dofPos[n]; it++) {
    if (it->type == dof) {
      pdof = &(*it);
      break;
    }
  }
  return pdof;
}

} // namespace nla3d

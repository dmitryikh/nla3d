#include "Dof.h"

namespace nla3d {

const uint16 Dof::numberOfDofTypes = Dof::UNDEFINED;
const char* const Dof::dofTypeLabels[] = {"UX", "UY", "UZ", "ROTX", "ROTY",
  "ROTZ", "HYDRO_PRESSURE", "UNDEFINED"};

uint32 DofCollection::empty = 0xFFFFFFFF;

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


void DofCollection::initDofTable(uint32 _numberOfEntities) {
  //TODO: make it possible to buildDofTable many times
  assert(dofPos.size() == 0);
  assert(dofs.size() == 0);

  numberOfEntities = _numberOfEntities;
  dofPos.assign(numberOfEntities + 1, empty);
  numberOfUsedDofs = 0;
}


// n index starts from 1
void DofCollection::addDof(uint32 n, Dof::dofType dof) {
  assert(n <= numberOfEntities);
  assert(dofPos.size() > 0);
  // if already registered - just return
  if (isDofUsed(n, dof)) return;

  if (dofPos[n] == empty) {
    // add to the end
    // if previus dofPos fields are empty we need to find last non empty and fill the gap
    uint32 last_no_zero = n-1;
    while (last_no_zero > 0 && dofPos[last_no_zero] == empty) last_no_zero--;
    uint32 lastPos = dofPos[last_no_zero];
    if (last_no_zero == 0 && lastPos == empty) lastPos = 0;
    assert(lastPos != empty);
    for (; last_no_zero < n; last_no_zero++) dofPos[last_no_zero] = lastPos;

    dofPos[n] = lastPos + 1;
    dofs.push_back(Dof(dof));
  } else {
    // need to insert dof.. bad case actualy
    // incement dofPos +1 for all above
    uint32 i = n;
    while (i <= numberOfEntities and dofPos[i] != empty) dofPos[i++]++;

    // insert Dof into dofs before dofPos[n]
    if (dofPos[n] == dofs.size()+1) {
      // can add to the end
      dofs.push_back(Dof(dof));
    } else if (dofPos[n] <= dofs.size()) {
      dofs.insert(dofs.begin() + (dofPos[n] - 1), Dof(dof));
    } else {
      assert(false);
    }

  }
  numberOfUsedDofs++;
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

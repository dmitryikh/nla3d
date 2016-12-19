// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"

namespace nla3d {

class DofCollection;
// class for Degree of Freedom informations 
// (is it fixed boundary condition, what its number in equation system, ..)
class Dof {
  public:
    enum dofType {
      UX = 0,
      UY,
      UZ,
      ROTX,
      ROTY,
      ROTZ,
      HYDRO_PRESSURE,
      UNDEFINED // should be last
    };

    uint32 eqNumber = 0; // 0 - not use
    bool isConstrained = false; // is this DOF defined by B.C.

    static const char* const dofTypeLabels[];
    static const uint16 numberOfDofTypes;
    static dofType label2dofType (const std::string& label);
    static const char* dofType2label (const dofType t);

    bool isUsed();

  friend class DofCollection;
  private:
    bool _isUsed = false;
};

inline const char* Dof::dofType2label(const Dof::dofType t) {
  assert(t < Dof::UNDEFINED);
  return dofTypeLabels[t];
}

inline bool Dof::isUsed() {
  return _isUsed;
}

class DofCollection {
  public:
    uint32 getNumberOfUsedDofs();
    uint32 getNumberOfAllocatedDofs();
    Dof* getDof(uint32 n, Dof::dofType dof);
    void buildDofTable(uint32 _numberOfEntities);
    void useDof(uint32 n, Dof::dofType dof);
    bool isDofUsed(uint32 n, Dof::dofType dof);
    void clearDofTable();

  private:
    uint32 numberOfUsedDofs = 0;
    uint32 numberOfAllocatedDofs = 0;
    uint32 numberOfEntities = 0;

    std::vector<Dof> dofTable;
};


inline uint32 DofCollection::getNumberOfUsedDofs() {
  return numberOfUsedDofs;
}


inline uint32 DofCollection::getNumberOfAllocatedDofs() {
  return numberOfAllocatedDofs;
}


// n index starts from 1
inline bool DofCollection::isDofUsed(uint32 n, Dof::dofType dof) {
  assert(n <= numberOfEntities);
  assert(dofTable.size() > 0);
  return dofTable[(n - 1) * Dof::numberOfDofTypes + dof]._isUsed;
}


inline Dof* DofCollection::getDof(uint32 n, Dof::dofType dof) {
  assert(n <= numberOfEntities);
  assert(dofTable.size() > 0);

  Dof* pdof = &(dofTable[(n - 1) * Dof::numberOfDofTypes + dof]);
  // return only used dofs
  assert(pdof->_isUsed == true);
  return pdof;
}

} // namespace nla3d 

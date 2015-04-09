// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"

namespace nla3d {

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
    UNDEFINED
  };

  static const char* const dofTypeLabels[];
  static const uint16 numberOfDofTypes;
  static dofType label2dofType (const std::string& label);

	Dof() : eqNumber(0), isConstrained(false) { }
	uint32 eqNumber; // 0 - not use
	bool isConstrained; // is this DOF defined by B.C.
};

} // namespace nla3d 

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once

namespace nla3d {
namespace query {

enum scalarQuery {
  SCALAR_UNDEF,
  SCALAR_SP,
  SCALAR_W,
  SCALAR_WU,
  SCALAR_WP,
  SCALAR_VOL,
  SCALAR_LAST
};

const char* const scalarQueryLabels[] = {"UNDEFINED", "S_P", "W", "WU", "WP", "VOL", "LAST"};


enum vectorQuery {
  VECTOR_UNDEF,
  VECTOR_IC,
  VECTOR_LAST
};

const char* const vectorQueryLabels[] = {"UNDEFINEDS", "IC", "LAST"};

enum tensorQuery {
	TENSOR_UNDEF,
  // usual stress tensor
	TENSOR_COUCHY, 
  // second Piola-Kirchgoff stress tensor (symmetric 3x3)
  TENSOR_PK2, 
  // Lagrange deformations
  TENSOR_E,  
  // C = F^T F 
  TENSOR_C,  
  TENSOR_LAST
};

const char* const tensorQueryLabels[] = {"UNDEFINED","COUCHY", "PK2", "E", "C"};

// means averaged values of the element
const uint16 GP_MEAN = 100;

} // namespace query
} // namespace nla3d

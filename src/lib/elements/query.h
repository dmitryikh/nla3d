// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once

namespace nla3d {

enum class scalarQuery {
  UNDEF = 0,
  SP,
  W,
  WU,
  WP,
  VOL,
  LAST
};

const char* const scalarQueryLabels[] = {"UNDEFINED", "S_P", "W", "WU", "WP", "VOL", "LAST"};

static_assert((int)scalarQuery::LAST == sizeof(scalarQueryLabels)/sizeof(scalarQueryLabels[0]) - 1,
    "scalarQuery enumeration and scalarQueryLabels must have the same number of entries");


enum class vectorQuery {
  UNDEF,
  IC,
  FLUX,
  GRADT,
  LAST
};

const char* const vectorQueryLabels[] = {"UNDEFINEDS", "IC", "FLUX", "GRADT", "LAST"};

static_assert((int)vectorQuery::LAST == sizeof(vectorQueryLabels)/sizeof(vectorQueryLabels[0]) - 1,
    "vectorQuery enumeration and vectorQueryLabels must have the same number of entries");


enum class tensorQuery {
	UNDEF,
  // usual stress tensor
	COUCHY, 
  // second Piola-Kirchgoff stress tensor (symmetric 3x3)
  PK2, 
  // Lagrange deformations
  E,  
  // C = F^T F 
  C,  
  LAST
};

const char* const tensorQueryLabels[] = {"UNDEFINED","COUCHY", "PK2", "E", "C", "LAST"};

static_assert((int)tensorQuery::LAST == sizeof(tensorQueryLabels)/sizeof(tensorQueryLabels[0]) - 1,
    "tensorQuery enumeration and tensorQueryLabels must have the same number of entries");

// means averaged values of the element
const uint16 GP_MEAN = 100;

inline char const* query2label(scalarQuery query) {
  assert(query >= scalarQuery::UNDEF && query < scalarQuery::LAST);
  return scalarQueryLabels[(int) query];
}

inline char const* query2label(vectorQuery query) {
  assert(query >= vectorQuery::UNDEF && query < vectorQuery::LAST);
  return vectorQueryLabels[(int) query];
}

inline char const* query2label(tensorQuery query) {
  assert(query >= tensorQuery::UNDEF && query < tensorQuery::LAST);
  return tensorQueryLabels[(int) query];
}

} // namespace nla3d

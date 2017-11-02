// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#pragma once
#include <string>
#include "math/Vec.h"
#include "math/Mat.h"

namespace nla3d {

class MaterialFactory;

// Material must have named material constants
class Material {
public:
  Material () : MC(NULL), numC(0), code(0) {
  }
  Material (uint16 num_c);
  ~Material()   // TODO: discover the virtual destructor
  {
    if (MC) delete[] MC;
    MC = NULL;
  }

  virtual std::string toString();

  std::string getName();
  uint16 getCode ();
  // constants getters
  double& Ci (uint16 i);
  double& Ci (const std::string& nameConst);
  uint16 getNumC ();

  static const double I[6];

  friend class MaterialFactory;
protected:
  void register_mat_const(uint16 num, ...);
  double* MC;
  std::vector<std::string> MC_names;
  uint16 numC;
  uint16 code;
  std::string name;
};


// ---=== FUNCTIONS ===--- //
inline std::string Material::getName()
{
  return name;
}

inline uint16 Material::getCode ()
{
  return code;
}

inline double& Material::Ci (uint16 i)
{
  assert(i < numC);
  return MC[i];
}

inline uint16 Material::getNumC ()
{
  return numC;
}

} // namespace nla3d

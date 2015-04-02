// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "materials/material.h"

namespace nla3d {

const double Material::I[6] = {1.0,0.0,0.0,1.0,0.0,1.0};

//---------------------------------------------------------
//----------------MATERIAL ABSTRACT CLASS------------------
//---------------------------------------------------------
Material::Material (uint16 num_c)
{
  code = 0;
	numC = num_c; 
	MC = new double[num_c];
}

std::string Material::toString() {
  std::string str;
	str += getName();
	for (size_t i=0; i < getNumC(); i++) {
	  str += " " + MC_names[i] + " = " + toStr(MC[i]);
  }
	return str;
}

double& Material::Ci (const std::string& nameConst) {
  for (size_t i = 0; i < getNumC(); i++) {
    if (nameConst.compare(MC_names[i]) == 0) {
      return MC[i];
    }
  }
	error("Material::Ci: can't find a material constant with name %s", nameConst.c_str());
}


void Material::register_mat_const(uint16 num, ...) {
	numC = num;
	va_list vlist;
	va_start(vlist, num);
  MC_names.clear();
  MC_names.reserve(numC);
	for (uint16 i=0; i < num; i++) {
		MC_names.push_back(va_arg(vlist,char*));
	}
	MC = new double[num];
}

} // namespace nla3d

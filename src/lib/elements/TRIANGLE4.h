// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"

namespace nla3d {

enum PlaneState{
  Strain,
  Stress
};

class ElementTRIANGLE4 : public ElementTRIANGLE {
public:
// ElementTRIANGLE4 () defines number of nodes in the element, number of dimensions (2D or 3D
// elements). It creates an array for storing global nodes numbers. And also, it registers which
// DoFs it is going to use (Dof::UX, Dof::UY, Dof::UZ in this case).
  ElementTRIANGLE4 ();

// pre() - functions that is called just before the solution procedures. pre() should register
// which element DoFs and nodal DoFs will be incorporated in the element stiffness matrix. On this
// step element alsoo need to initialize any variables that it is going to use in solution process
// (strains and stresses in integration points in finite deformations analysis, for example).
// ElementTRIANGLE4::pre () registers Dof::UX, Dof::UY, Dof::UZ as DoFs in every node.
  void pre();

// buildK() - a central point in element class. Here the element should build element stiffness
// matrix (actually, tangential matrix, as soon as we make non-linear-ready elements). The element
// also responsible for assembling its local stiffness matrix into global system of equations
// matrix. Fro this purpose here is a special procedure in base class Element::assemble(..). Also,
// the element should assemble right hand side (rhs) of equations related to this element
// (especially used in non-linear analysis).
//
  void buildK();
// update() - the function updates internal state of the element based on found solution of
// global equation system. For example, here you can calculate stresses in the element which depends
// on found DoFs solution.
  void update();


  void makeB (math::Mat<3,6> &B);
  void makeC (math::MatSym<3> &C);

  // Elastic module
  double E = 0.0;
  // Poissons coef.
  double my = 0.0;

  // stresses in the element (calculated after the solving of the global equation system in
  // update() function.
  //stress[M_X], stress[M_X], stress[M_XY]
  math::Vec<3> stress; // Cauchy stresses

  //strains[M_X], strains[M_Y], strains[M_XY]
  math::Vec<3> strains;

  PlaneState state;
private:
  double area;
};

} //namespace nla3d

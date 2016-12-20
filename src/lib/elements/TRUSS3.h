// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"

namespace nla3d {

// ElementTRUSS3 (TRUSS3) is a simplest realisation of 3D truss finite element (FE). Its primiral goal is to
// show how to implement new FE into nla3d library. One can find a lot of information about how to
// build 3D truss element stiffnes matrix in the internet. The author was looking at the material on
// this website: 
// http://what-when-how.com/the-finite-element-method/fem-for-trusses-finite-element-method-part-1/

// New FE formulation is encapsulated into a class derived from base Element class. FEs have a lot
// of freedom in nla3d. They are responsible for registration degrees of freedom (DoFs) that they
// will use (see pre() funciton). They fully responsible how to work with materials instances. The
// class ElelementTRUSS3 doesn't use material information at all. Instead of this it uses class
// members A and E as truss properties. The author thinks that FE and material can't be fully
// divided into different parts, as soon as FE stiffness matrix formulation highly depends on which
// material we use (elastic, anisotropic, plastic, hyperelastic, etc). Then elements in nla3d have a
// right to not use Material class (see src/lib/materials/material.h) at all.

class ElementTRUSS3 : public ElementLINE {
public:
// Element in nla3d have to provide several obligatore functions in order to make it possible to use
// it:
// 1. element class constructor to define number of nodes, dimension, and register which DoF types the element
// are going to use.
// ElementTRUSS3 () defines number of nodes in the element, number of dimensions (2D or 3D
// elements). It creates an array for storing global nodes numbers. And also, it registers which
// DoFs it is going to use (Dof::UX, Dof::UY, Dof::UZ in this case).
  ElementTRUSS3 ();
// 2. pre() - functions that is called just before the solution procedures. pre() should register
// which element DoFs and nodal DoFs will be incorporated in the element stiffness matrix. On this
// step element alsoo need to initialize any variables that it is going to use in solution process
// (strains and stresses in integration points in finite deformations analysis, for example).
// ElementTRUSS3::pre () registers Dof::UX, Dof::UY, Dof::UZ as DoFs in every node.
  void pre();
// 3. build() - a central point in element class. Here the element should build element stiffness
// matrix (actually, tangential matrix, as soon as we make non-linear-ready elements). The element
// also responsible for assembling its local stiffness matrix into global system of equations
// matrix. Fro this purpose here is a special procedure in base class Element::assemble(..). Also,
// the element should assemble right hand side (rhs) of equations related to this element
// (especially used in non-linear analysis).
// see ElementTRUSS3::build() body for more comments on the particular realisation.
  void build();
// 4. update() - the function updates internal state of the element based on found solution of
// global equation system. For example, here you can calculate stresses in the element which depends
// on found DoFs solution.
// see ElementTRUSS3::update() body for the insight of the particular realization.
  void update();
  // stiffness module
  double E = 0.0;
  // cross-section area
  double A = 0.0;
  // normal stress in the truss (calculated after the solving of the global equation system in
  // update() function.
  double S;
};

} //namespace nla3d

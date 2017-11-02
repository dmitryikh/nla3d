// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "sys.h"
#include "FEStorage.h"
#include "VtkProcessor.h"
#include "FESolver.h"
#include "elements/TRUSS3.h"

using namespace nla3d;

// This exceutable program is a continuation of educational TRUSS3 element realization. As it was
// said, nla3d is not a ready-to-use FE program. It's a library of routines to implement new FE
// algorithms. After ElementTRUSS3 FE was written, one needs to have an exceutable to solve some
// FE model using this element. Here is an example of exceutable to solve a particular truss system.
// The description of this truss system you can found in files/truss-method-example.pdf (on page 4).
// Originally this file was downloaded from http://people.duke.edu/~hpgavin/cee421/truss-method.pdf.
// In the code below you will find all neede comments about how all this work.

int main (int argc, char* argv[]) {
  // Definition of the node table. Every node in nla3d lives in 3D space (has 3 coordinates). As fas
  // as we have a deal with 2D case, we just leave third coordinate equal to zero.
  const uint32 numberOfNodes = 5;
  double nodeTable[numberOfNodes][3] = {
                          {0.0, 0.0, 0.0},
                          {16.0, 0.0, 0.0},
                          {16.0, 12.0, 0.0},
                          {32.0, 0.0, 0.0},
                          {32.0, 12.0, 0.0}};

  // Definition of the element table. It's clear that TRUSS3 element has 2 nodes.
  const uint32 numberOfElements = 8;
  uint32 elementTable[numberOfElements][2] = {
                              {1, 3},
                              {1, 2},
                              {2, 3},
                              {3, 5},
                              {3, 4},
                              {2, 5},
                              {2, 4},
                              {4, 5}};

  // Create an instance of FEStorage.
	FEStorage storage;

  // We have a deal with linear FE. Then it's ok to use linear solver (just one equilibrium iteration without
  // convergence controls)
	LinearFESolver solver;

  // Create and add nodes into FEStorage
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    Node* no = new Node;
    no->pos[0] = nodeTable[i-1][0]*12.0;
    no->pos[1] = nodeTable[i-1][1]*12.0;
    no->pos[2] = nodeTable[i-1][2]*12.0;
    storage.addNode(no);
  }

  // Create elements instances, define needed element parameters and add them into FEStorage.
  for (uint32 i = 1; i <= numberOfElements; i++) {
    ElementTRUSS3* el = new ElementTRUSS3;
    el->E = 3.0e4;
    el->A = 10.0;
    el->getNodeNumber(0) = elementTable[i-1][0];
    el->getNodeNumber(1) = elementTable[i-1][1];
    storage.addElement(el);
  }

  // we deal with 2D system, then fix all UZ DoFs
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    solver.addFix(i, Dof::UZ);
  }

  // add fixations for the system
  solver.addFix(1, Dof::UX);
  solver.addFix(1, Dof::UY);
  solver.addFix(4, Dof::UX);
  solver.addFix(4, Dof::UY);

  // add forces
  solver.addLoad(2, Dof::UY, -100.0);
  solver.addLoad(5, Dof::UX, 50.0);

#ifdef NLA3D_USE_MKL
    math::PARDISO_equationSolver eqSolver = math::PARDISO_equationSolver();
    solver.attachEquationSolver(&eqSolver);
#endif
  // FESolver should know FEStorage instance. Attach it.
	solver.attachFEStorage(&storage);
  // We would like to generate *.vtk file with deformed and undeformed models. For this purpose
  // PostProcessor with name VtkProcessor is added to FESolver. The names of vtk files will be
  // "truss2D0.vtk" for undeformed model, and "truss2D1.vtk" for deformed one.
  auto vtk = std::make_shared<VtkProcessor>(&storage, "truss2D");
  solver.addPostProcessor(vtk);
  // just solve the model. Yes, so easy.
	solver.solve();
  
  // Log all results about the model
  LOG(INFO) << "DoF solution:";
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    LOG(INFO) << i << ":" << Dof::dofTypeLabels[Dof::UX] << " = " << storage.getNodeDofSolution(i, Dof::UX);
    LOG(INFO) << i << ":" << Dof::dofTypeLabels[Dof::UY] << " = " << storage.getNodeDofSolution(i, Dof::UY);
  }

  LOG(INFO) << "DoF reactions:";
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    LOG(INFO) << i << ":" << Dof::dofTypeLabels[Dof::UX] << " = " << storage.getReaction(i, Dof::UX);
    LOG(INFO) << i << ":" << Dof::dofTypeLabels[Dof::UY] << " = " << storage.getReaction(i, Dof::UY);
  }

  LOG(INFO) << "Stresses in elements:";
  for (uint32 i = 1; i <= numberOfElements; i++) {
    ElementTRUSS3* el = dynamic_cast<ElementTRUSS3*> (&storage.getElement(i));
    LOG(INFO) << i << " = " << el->S;
  }
  // check for correct result
  CHECK(storage.getNodeDofSolution(2, Dof::UX) - 0.0146067 < 1.0e-7);
  CHECK(storage.getNodeDofSolution(2, Dof::UY) - (-0.1046405) < 1.0e-7);

  CHECK(storage.getNodeDofSolution(3, Dof::UX) - 0.0027214 < 1.0e-7);
  CHECK(storage.getNodeDofSolution(3, Dof::UY) - (-0.0730729) < 1.0e-7);

  CHECK(storage.getNodeDofSolution(5, Dof::UX) - 0.0055080 < 1.0e-7);
  CHECK(storage.getNodeDofSolution(5, Dof::UY) - (-0.0164325) < 1.0e-7);
	return 0;
}

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "sys.h"
#include "FEStorage.h"
#include "VtkProcessor.h"
#include "FESolver.h"
#include "elements/QUADTH.h"

using namespace nla3d;


std::vector<double> readTempData (std::string fileRefCurve);

int main (int argc, char* argv[]) {

  std::string cdb_filename = "";
  std::string res_filename = "";

  if (argc > 1) {
    cdb_filename = argv[1];
  } else {
    LOG(FATAL) << "You should provide the path to mesh (cdb file)";
  }

  if (argc > 2) {
    res_filename = argv[2];
  }

  // Definition of the element table for SurfaceLINETH
  const uint32 numberOfElements = 5;
  uint32 elementTable[numberOfElements][2] = {
                              {6, 5},
                              {5, 4},
                              {4, 3},
                              {3, 1},
                              {2, 6}};

  // Create an instance of FEStorage.
	FEStorage storage;
  if (!readCdbFile (cdb_filename.c_str(), &storage, ElementType::QUADTH)) {
    LOG(ERROR) << "Can't read FE info from cdb";
    exit(1);
  }

  for (uint32 i = 1; i <= storage.getNumberOfElements(); i++) {
    ElementQUADTH& el = dynamic_cast<ElementQUADTH&>(storage.getElement(i));
    el.k = 0.018;
  }

  // Create elements instances, define needed element parameters and add them into FEStorage.
  for (uint32 i = 1; i <= numberOfElements; i++) {
    SurfaceLINETH* el = new SurfaceLINETH;
    el->htc = 0.0034;
    el->etemp[0] = -6.0;
    el->etemp[1] = -6.0;
    el->getNodeNumber(0) = elementTable[i-1][0];
    el->getNodeNumber(1) = elementTable[i-1][1];
    storage.addElement(el);
  }

  storage.addDofLoad(7, Dof::TEMP, 0.08);

  // We have a deal with linear FE. Then it's ok to use linear solver (just one equilibrium iteration without
  // convergence controls)
	LinearFESolver solver;
#ifdef NLA3D_USE_MKL
    math::PARDISO_equationSolver eqSolver = math::PARDISO_equationSolver();
    solver.attachEquationSolver(&eqSolver);
#endif
  // FESolver should know FEStorage instance. Attach it.
	solver.attachFEStorage(&storage);
  VtkProcessor* vtk = new VtkProcessor(&storage, "QUADTH");
  solver.addPostProcessor(vtk);

	solver.solve();
  
  // Log all results about the model
  LOG(INFO) << "DoF solution:";
  for (uint32 i = 1; i <= storage.getNumberOfNodes(); i++) {
    LOG(INFO) << i << ":" << Dof::dofTypeLabels[Dof::TEMP] << " = " << storage.getNodeDofSolution(i, Dof::TEMP);
  }


  // check results with Ansys data
  if (res_filename != "") {
    auto ansTemp = readTempData (res_filename);
    for (uint32 i = 1; i <= storage.getNumberOfNodes(); i++) {
      double my_temp = storage.getNodeDofSolution(i, Dof::TEMP);
      double ans_temp = ansTemp[i-1];
      CHECK_EQTH(my_temp, ans_temp, 1.0e-3);
    }
  }

	return 0;
}


std::vector<double> readTempData (std::string fileRefCurve) {
  std::ifstream file(fileRefCurve.c_str());
  double tmp;
  std::vector<double> res;
  // skip first line
  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  while (!file.eof()) {
    // read second column
    file >> tmp >> tmp;
    res.push_back(tmp);
  }
  file.close();
  return res;
}

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "sys.h"
#include <limits>
#include "FEStorage.h"
#include "VtkProcessor.h"
#include "FESolver.h"
#include "elements/QUADTH.h"
#include "FEReaders.h"

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

  MeshData md;
  if (!readCdbFile(cdb_filename, md)) {
    LOG(FATAL) << "Can't read FE data from cdb";
  }
  md.compressNumbers();

  if (md.getDegeneratedCells().size() > 0) {
    LOG(FATAL) << "Can't work with degenerated element in the mesh";
  }

  // Create an instance of FEStorage.
	FEStorage storage;
  // We have a deal with linear FE. Then it's ok to use linear solver (just one equilibrium iteration without
  // convergence controls)
	LinearFESolver solver;

  auto sind = storage.createNodes(md.nodesNumbers.size());
  auto ind = md.nodesNumbers;
  for (uint32 i = 0; i < sind.size(); i++) {
    storage.getNode(sind[i]).pos = md.nodesPos[i];
  }

  ind = md.getCellsByAttribute("TYPE", 1);
  sind = storage.createElements(ind.size(), ElementQUADTH());
  for (uint32 i = 0; i < sind.size(); i++) {
    ElementQUADTH& el = dynamic_cast<ElementQUADTH&>(storage.getElement(sind[i]));
    assert(el.getNNodes() == md.cellNodes[ind[i]].size());
    for (uint16 j = 0; j < el.getNNodes(); j++) {
      el.getNodeNumber(j) = md.cellNodes[ind[i]][j];
    }
    el.k = 0.018;
  }


  ind = md.getCellsByAttribute("TYPE", 2);
  sind = storage.createElements(ind.size(), SurfaceLINETH());
  for (uint32 i = 0; i < sind.size(); i++) {
    SurfaceLINETH& el = dynamic_cast<SurfaceLINETH&>(storage.getElement(sind[i]));
    assert(el.getNNodes() == md.cellNodes[ind[i]].size());
    for (uint16 j = 0; j < el.getNNodes(); j++) {
      el.getNodeNumber(j) = md.cellNodes[ind[i]][j];
    }
    el.htc = 0.0034;
    el.etemp[0] = -6.0;
    el.etemp[1] = -6.0;
  }

  // add loadBc
  for (auto& v : md.loadBcs) {
    solver.addLoad(v.node, v.node_dof, v.value);
  }

  // add fixBc
  for (auto& v : md.fixBcs) {
    solver.addFix(v.node, v.node_dof, v.value);
  }

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
  for (uint32 i = 1; i <= storage.nNodes(); i++) {
    LOG(INFO) << i << ":" << Dof::dofTypeLabels[Dof::TEMP] << " = " << storage.getNodeDofSolution(i, Dof::TEMP);
  }


  // check results with Ansys data
  if (res_filename != "") {
    auto ansTemp = readTempData (res_filename);
    for (uint32 i = 1; i <= storage.nNodes(); i++) {
      double my_temp = storage.getNodeDofSolution(i, Dof::TEMP);
      double ans_temp = ansTemp[i-1];
      CHECK_EQTH(my_temp, ans_temp, 1.0e-3);
    }
  }

	return 0;
}


std::vector<double> readTempData (std::string fileRefCurve) {
  std::vector<double> res;
  std::ifstream file(fileRefCurve.c_str());
  if (!file.is_open()) {
    return res;
  }
  double tmp;
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

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "sys.h"
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
	LinearTransientFESolver solver;

  auto sind = storage.createNodes(md.nodesNumbers.size());
  auto ind = md.nodesNumbers;
  for (uint32 i = 0; i < sind.size(); i++) {
    storage.getNode(sind[i]).pos = md.nodesPos[i];
  }

  // take elements with TYPE = 1 from MeshData and create corresponding ElementQUADTH elements in
  // FEStorage
  ind = md.getCellsByAttribute("TYPE", 1);
  sind = storage.createElements(ind.size(), ElementQUADTH());
  for (uint32 i = 0; i < sind.size(); i++) {
    ElementQUADTH& el = dynamic_cast<ElementQUADTH&>(storage.getElement(sind[i]));
    assert(el.getNNodes() == md.cellNodes[ind[i]].size());
    for (uint16 j = 0; j < el.getNNodes(); j++) {
      el.getNodeNumber(j) = md.cellNodes[ind[i]][j];
    }
    // setup element constants
    el.k = 0.018; // W/cm C
    el.c = 920; // J/kg C
    el.rho = 1.1e-3; // kg/cm³
  }

  // TYPE = 2 elements it's Surface (line in 2D case) elements to implement convection BC
  ind = md.getCellsByAttribute("TYPE", 2);
  sind = storage.createElements(ind.size(), SurfaceLINETH());
  for (uint32 i = 0; i < sind.size(); i++) {
    SurfaceLINETH& el = dynamic_cast<SurfaceLINETH&>(storage.getElement(sind[i]));
    assert(el.getNNodes() == md.cellNodes[ind[i]].size());
    for (uint16 j = 0; j < el.getNNodes(); j++) {
      el.getNodeNumber(j) = md.cellNodes[ind[i]][j];
    }
    el.htc = 0.0034; // wall heat transfer coefficient
    el.etemp[0] = -6.0; // bulk temperature at 1-st node
    el.etemp[1] = -6.0; // bulk temperature at 2-nd node
  }

  // component C0 contains single node for applying thermal source
  solver.addLoad(md.feComps["C0"].list[0], Dof::TEMP, 0.08);

  // time period to solve in seconds
  solver.time1 = 30000.0;
  solver.numberOfTimesteps = 100;
  // value to initialize initial field
  solver.initValue = -6.0;
#ifdef NLA3D_USE_MKL
    math::PARDISO_equationSolver eqSolver = math::PARDISO_equationSolver();
    solver.attachEquationSolver(&eqSolver);
#endif
  // FESolver should know FEStorage instance. Attach it.
	solver.attachFEStorage(&storage);
  VtkProcessor* vtk = new VtkProcessor(&storage, "QUADTH_transient");
  solver.addPostProcessor(vtk);

	solver.solve();
  
  // Log all results about the model
  LOG(INFO) << "DoF solution:";
  for (uint32 i = 1; i <= storage.nNodes(); i++) {
    LOG(INFO) << i << ":" << Dof::dofTypeLabels[Dof::TEMP] << " = " << storage.getNodeDofSolution(i, Dof::TEMP);
  }


  // check results with Ansys data
  // if (res_filename != "") {
  //   auto ansTemp = readTempData (res_filename);
  //   for (uint32 i = 1; i <= storage.getNumberOfNodes(); i++) {
  //     double my_temp = storage.getNodeDofSolution(i, Dof::TEMP);
  //     double ans_temp = ansTemp[i-1];
  //     CHECK_EQTH(my_temp, ans_temp, 1.0e-3);
  //   }
  // }

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

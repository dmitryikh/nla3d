// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d
//
#include "sys.h"
#include "FEStorage.h"
#include "FESolver.h"
#include "VtkProcessor.h"
#include "ReactionProcessor.h"
#include "materials/MaterialFactory.h"
#include "FEReaders.h"
#include "elements/TETRA1.h"
#include <tuple>

using namespace nla3d;

typedef std::vector<double> temp_vec_t;

temp_vec_t readTermData (std::string filename);

int main (int argc, char* argv[]) {
  std::string cdb_filename;
  std::string res_term_filename;

  if (argc > 1) {
    cdb_filename = argv[1];
  } else {
    LOG(FATAL) << "You should provide the path to mesh (cdb file)";
  }
  if (argc > 2) {
    res_term_filename = argv[2];
  }

  MeshData md;
  if (!readCdbFile(cdb_filename, md)) {
    LOG(FATAL) << "Can't read FE info from " << cdb_filename << "file. exiting..";
  }
  md.compressNumbers();

  FEStorage storage;
  LinearFESolver solver;

  // add nodes
  auto sind = storage.createNodes(md.nodesNumbers.size());
  auto ind = md.nodesNumbers;
  for (uint32 i = 0; i < sind.size(); i++) {
    storage.getNode(sind[i]).pos = md.nodesPos[i];
  }

  ind = md.getCellsByAttribute("TYPE", 1);
  sind = storage.createElements(ind.size(), ElementType::TETRA1);
  for (uint32 i = 0; i < sind.size(); i++) {
    ElementTETRA1& el = dynamic_cast<ElementTETRA1&>(storage.getElement(sind[i]));
    el.getNodeNumber(0) = md.cellNodes[ind[i]][0];
    el.getNodeNumber(1) = md.cellNodes[ind[i]][1];
    el.getNodeNumber(2) = md.cellNodes[ind[i]][2];
    el.getNodeNumber(3) = md.cellNodes[ind[i]][4];
    el.k = 60.5; //W*m-1*C-1
  }

  for (auto v : md.feComps["X0"].list) {
    solver.addFix(v, Dof::TEMP,0.);
  }
  for (auto v : md.feComps["X1"].list) {
    solver.addFix(v, Dof::TEMP,10.);
  }
  for (auto v : md.feComps["Y"].list) {
    solver.addFix(v, Dof::TEMP,5.);
  }
  for (auto v : md.feComps["ORIGIN"].list) {
    solver.addFix(v, Dof::TEMP,20.);
  }

#ifdef NLA3D_USE_MKL
  math::PARDISO_equationSolver eqSolver = math::PARDISO_equationSolver();
  solver.attachEquationSolver(&eqSolver);
#endif
  solver.attachFEStorage (&storage);

  VtkProcessor* vtk = new VtkProcessor (&storage, "tetra1");
  solver.addPostProcessor(vtk);
  vtk->writeAllResults();

  solver.solve();

  // check terperature results with Ansys data
  if (res_term_filename != "") {
    auto ans_temp_vec = readTermData(res_term_filename);
    for (uint32 i = 1; i <= storage.nNodes(); i++) {
      const double t = storage.getNodeDofSolution(i, Dof::TEMP);
      const auto ans_temp = ans_temp_vec[i - 1];
      CHECK_EQTH(t, ans_temp, 1.0e-3);
    }
  }
}

temp_vec_t readTermData (std::string filename) {
  temp_vec_t tepm_vec;
  std::ifstream file(filename);
  if (!file.is_open()) {
    return tepm_vec;
  }
  double node, temp;
  // file format:
  //Node Number Temperature (°C)
  //1   20.
  //  ...

  // skip first line
  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  while (!file.eof()) {
    file >> node >> temp;
    tepm_vec.emplace_back(temp);
  }
  return tepm_vec;
}
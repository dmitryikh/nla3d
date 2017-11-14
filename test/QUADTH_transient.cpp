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

class ProbeProcessor : public PostProcessor {
public:
	ProbeProcessor(FEStorage *st);
	virtual ~ProbeProcessor() { };

	virtual void pre ();
	virtual void process(uint16 curLoadstep);
	virtual void post(uint16 curLoadstep);

  std::vector<double> getMeasurments();
	
	uint32 node = 0;
	Dof::dofType dof = Dof::UNDEFINED;
protected:
	std::vector<double> measurments;
};


ProbeProcessor::ProbeProcessor(FEStorage *st) : PostProcessor(st) {
	name ="ProbeProcessor";
}


void ProbeProcessor::pre() {
	if (node == 0 || dof == Dof::UNDEFINED) {
    LOG(FATAL) << "Can't work. No node number and/or Dof id. Processor name = " << name;
    return;
  }
}


void ProbeProcessor::process(uint16 curLoadstep) {
  double val = storage->getNodeDofSolution(node, dof);
  measurments.push_back(val);
}


void ProbeProcessor::post(uint16 curLoadstep) {
}


std::vector<double> ProbeProcessor::getMeasurments() {
  return measurments;
}

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

  // add loadBc
  for (auto& v : md.loadBcs) {
    solver.addLoad(v.node, v.node_dof, v.value);
  }

  // add fixBc
  for (auto& v : md.fixBcs) {
    solver.addFix(v.node, v.node_dof, v.value);
  }

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

  ProbeProcessor* probe = new ProbeProcessor(&storage);
  probe->node = md.feComps["PROBE"].list[0];
  probe->dof = Dof::TEMP;
  solver.addPostProcessor(probe);

	solver.solve();
  
  // Log all results about the model
  LOG(INFO) << " Probe Temperature:";
  auto meas = probe->getMeasurments();
  for (uint32 i = 1; i <= meas.size(); i++) {
    LOG(INFO) << float(i * solver.time1) / float(solver.numberOfTimesteps)  << ": " << meas[i - 1];
  }


  // check results with Ansys data
  if (res_filename != "") {
    double ave_fabs = 0.0;
    auto ansTemp = readTempData(res_filename);
    for (uint32 i = 1; i <= meas.size(); i++) {
      double my_temp = meas[i - 1];
      double ans_temp = ansTemp[i - 1];
      ave_fabs += fabs(my_temp - ans_temp);
    }
    ave_fabs /= meas.size();
    LOG(INFO) << "Average error is: " << ave_fabs;
    CHECK(ave_fabs < 5.0e-2);
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
  std::string line;
  // skip first line
  getLine(file, line);

  Tokenizer t;
  t.delimiters.insert(',');

  while (!getLine(file, line).eof()) {
    // read second column
    t.tokenize(line);
    res.push_back(t.tokenDouble(1));
  }
  file.close();
  return res;
}

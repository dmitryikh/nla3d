// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "sys.h"
#include "FEStorage.h"
#include "VtkProcessor.h"
#include "Solution.h"
#include "Reaction_proc.h"
#include "materials/MaterialFactory.h"

using namespace nla3d;

void prepeareProcessors (FEStorage& storage);
std::vector<double> readRefCurveData (std::string fileRefCurve); 
double compareCurves (const std::vector<double>& refCurve, const std::vector<double>& curCurve);

// Here is defaults values for command line options
// For command line format see usage() below
namespace options {
  uint16 numberOfIterations = 15;
  uint16 numberOfLoadsteps = 10;
  std::string materialName = "";
  ElementFactory::elTypes elementType = ElementFactory::SOLID81;
  bool useVtk = true;
  std::string modelFilename = "";
  std::vector<double> materialConstants;
  std::string refCurveFilename = ""; 
  double curveCompareThreshold = 0.0;
  std::string reactionComponentName = "";
};

bool parse_args (int argc, char* argv[]) {
  if (argc < 2) {
    return false;
  }
  options::modelFilename = argv[1];
  argv = argv + 2;
  argc = argc - 2;
  if(cmdOptionExists(argv, argv+argc, "-novtk")) {
    options::useVtk = false;
  }

  char* tmp = getCmdOption(argv, argv + argc, "-iterations");
  if (tmp) {
    options::numberOfIterations = atoi(tmp);
  }

  tmp = getCmdOption(argv, argv + argc, "-loadsteps");
  if (tmp) {
    options::numberOfLoadsteps = atoi(tmp);
  }

  tmp = getCmdOption(argv, argv + argc, "-element");
  if (tmp) {
    options::elementType = ElementFactory::elName2elType(tmp);
  }

  std::vector<char*> vtmp = getCmdManyOptions(argv, argv + argc, "-material");
  if (vtmp.size() == 0) {
    error("Please point a material model. Use -material keyword.");
  } else {
    options::materialName = vtmp[0];
    for (uint16 i = 1; i < vtmp.size(); i++) {
      options::materialConstants.push_back(atof(vtmp[i]));
    }
  }

  tmp = getCmdOption(argv, argv + argc, "-refcurve");
  if (tmp) {
    options::refCurveFilename = tmp;
  }

  tmp = getCmdOption(argv, argv + argc, "-threshold");
  if (tmp) {
    options::curveCompareThreshold = atof(tmp);
  }

  tmp = getCmdOption(argv, argv + argc, "-reaction");
  if (tmp) {
    options::reactionComponentName = tmp;
  }

  return true;
}

void usage () {
  echolog("nla3d model_file [-element el_name] [-material mat_name mat_C0 mat_C1 ..] [-iterations it_num] [-loadsteps ls_num] [-novtk]");
}

int main (int argc, char* argv[])
{

  std::vector<double> refCurve;
  std::vector<double> curCurve;

  Reaction_proc* reactProc;
	echolog("---=== WELCOME TO NLA PROGRAM ===---");
  if (!parse_args(argc, argv)) {
    usage();
    exit(1);
  }

  // Load and show reference curve. This is used for automatic tests
  // and optimization (curve fitting) purpose.
  if (options::refCurveFilename.length() > 0) {
    if (options::reactionComponentName.length() == 0) {
      error("Reference curve is obtained, but no component name provided to \
          calculate reactions from analysis. Use -reaction option");
    }
    refCurve = readRefCurveData(options::refCurveFilename);
    for (size_t i = 0; i < refCurve.size(); i++) {
      echolog("refCurve[%d] = %f", i, refCurve[i]);
    }
  }

	Timer pre_solve(true);
	FEStorage storage;
  storage.elType = options::elementType;
  if (!read_ans_data(options::modelFilename.c_str(), &storage)) {
    error("Can't read FE info from %s file. exiting..", options::modelFilename.c_str());
  }
  Material* mat = MaterialFactory::createMaterial(options::materialName);
  if (mat->getNumC() != options::materialConstants.size())
    error("Material %s needs exactly %d constants (%d were provided)", options::materialName.c_str(),
        mat->getNumC(), options::materialConstants.size());
  for (uint16 i = 0; i < mat->getNumC(); i++) {
    mat->Ci(i) = options::materialConstants[i];
  }
	storage.material = mat;
  echolog("Material: %s", mat->toString().c_str());

  std::stringstream ss;
  ss << "Loaded components:" << std::endl;
  storage.listFEComponents(ss);
  echolog(ss.str().c_str());

	Solution sol;
	sol.attach(&storage);
	sol.setqIterat(options::numberOfIterations);
	sol.setqLoadstep(options::numberOfLoadsteps);
  if (options::useVtk) {
    //obtain job name from path of a FE model file
    std::string jobname = getFileNameFromPath(options::modelFilename);
    VtkProcessor* vtk = new VtkProcessor (&storage, jobname);
  }

  reactProc = NULL;
  if (options::reactionComponentName.length() > 0) {
    reactProc = new Reaction_proc(&storage);
    FEComponent* feComp = storage.getFEComponent(options::reactionComponentName);
    if (feComp->type != FEComponent::NODES) {
      error("To calculate reactions it's needed to provide a component with nodes");
    }

    echolog("Component for measuring a reaction: %s", options::reactionComponentName.c_str());
    for (size_t i = 0; i < feComp->list.size(); i++) {
      reactProc->nodes.push_back(feComp->list[i]);
      reactProc->dofs.push_back(1);
    }
  }

	echolog("Preprocessor time: %f sec.", pre_solve.stop());
	sol.run();

  if (reactProc) {
    curCurve = reactProc->getReactions(); 
    for (size_t i = 0; i < curCurve.size(); i++) {
      echolog("curCurve[%d] = %f", i, curCurve[i]);
    }
    if (refCurve.size() > 0) {
      double _error = compareCurves (refCurve, curCurve);
      echolog("Error between reference loading curve and current is %f", _error);
      if (options::curveCompareThreshold > 0.0) {
        if (_error > options::curveCompareThreshold) {
          error("To big error! (upper bound is %f)", options::curveCompareThreshold);
        } else {
          echolog("Error is less that the threshold. %f < %f", _error, options::curveCompareThreshold);
        }
      }
    }
  }

	return 0;
}


double compareCurves (const std::vector<double>& refCurve, const std::vector<double>& curCurve) {
  if (refCurve.size() != curCurve.size()) {
    error("curves: tabular data size should be the same for reference curve \
        and for current curve. (got refCurve.size() = %d, curCurve.size() = %d",
        refCurve.size(), curCurve.size());
  }
  double _error = 0.0;
  size_t lenCurve = curCurve.size(); 
  for (size_t i = 0; i < lenCurve; i++) {
    _error += fabs(refCurve[i]-curCurve[i])/lenCurve;
  }
  return _error;
}

std::vector<double> readRefCurveData (std::string fileRefCurve) {
  std::ifstream file(fileRefCurve.c_str());
  std::string dummy;
  double tmp;
  double pre_tmp = -1.0;
  std::vector<double> res;
  file >> dummy >> dummy >> dummy;
  if (dummy.compare("force") != 0) {
    error("readRefCurveData: Expecting third column to a force column");
  }
  while (!file.eof()) {
    file >> dummy >> dummy >> tmp;
    if (tmp != pre_tmp) {
      res.push_back(tmp);
      pre_tmp = tmp;
    }
  }
  file.close();
  return res;
}

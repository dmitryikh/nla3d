// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "sys.h"
#include "FEStorage.h"
#include "FESolver.h"
#include "VtkProcessor.h"
#include "ReactionProcessor.h"
#include "materials/MaterialFactory.h"
#include "FEReaders.h"

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
  ElementType elementType = ElementType::SOLID81;
  bool useVtk = true;
  std::string modelFilename = "";
  std::vector<double> materialConstants;
  std::string refCurveFilename = ""; 
  double curveCompareThreshold = 0.0;
  std::string reactionComponentName = "";
  std::vector<Dof::dofType> reactionDofs;

  uint32 rigidBodyMasterNode = 0;
  std::string rigidBodySlavesComponent = "";
  std::vector<Dof::dofType> rigidBodyDofs;
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
    LOG(ERROR) << "Please point a material model. Use -material keyword.";
    std::exit(1);
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

  vtmp = getCmdManyOptions(argv, argv + argc, "-reaction");
  if (vtmp.size() > 0) {
    if (vtmp.size() == 1) {
      options::reactionComponentName = vtmp[0];
    } else  {
      options::reactionComponentName = vtmp[0];
      for (size_t i = 1; i < vtmp.size(); i++) {
        options::reactionDofs.push_back(Dof::label2dofType(vtmp[i])); 
      }
    }
  }

  vtmp = getCmdManyOptions(argv, argv + argc, "-rigidbody");
  if (vtmp.size() > 0) {
    if (vtmp.size() < 2) {
      LOG(ERROR) << "Please point master node, componen with slave node and degrees of freedom";
      std::exit(1);
    } else {
      options::rigidBodyMasterNode = atoi(vtmp[0]);
      options::rigidBodySlavesComponent = vtmp[1];
      if (vtmp.size() == 2) {
        options::rigidBodyDofs.push_back(Dof::UX);
        options::rigidBodyDofs.push_back(Dof::UY);
        options::rigidBodyDofs.push_back(Dof::UZ);
      } else {
        for (size_t i = 2; i < vtmp.size(); i++) {
          options::rigidBodyDofs.push_back(Dof::label2dofType(vtmp[i])); 
        }
      }
    }
  }

  return true;
}

void usage () {
  LOG(INFO) << "nla3d 'mesh file in cdb format'\n"
      << "\t[-element 'element name']\n"
      << "\t[-material 'material name' constant1 constant2 ..]\n"
      << "\t[-iterations 'number of iterations']\n"
      << "\t[-loadsteps 'number of loadsteps']\n"
      << "\t[-novtk]\n"
      << "\t[-refcurve 'file with curve']\n"
      << "\t[-threshold 'epsilob for comparison']\n"
      << "\t[-reaction 'component name' ['DoF' ..]]\n"
      << "\t[-rigidbody 'master node' 'component of slaves' ['DoF' ..]]";
}

int main (int argc, char* argv[]) {
  std::vector<double> refCurve;
  std::vector<double> curCurve;

  ReactionProcessor* reactProc;
  LOG(INFO) << "---=== WELCOME TO NLA PROGRAM ===---";
  if (!parse_args(argc, argv)) {
    usage();
    exit(1);
  }

  // Load and show reference curve. This is used for automatic tests
  // and optimization (curve fitting) purpose.
  if (options::refCurveFilename.length() > 0) {
    if (options::reactionComponentName.length() == 0) {
      LOG(ERROR) << "Reference curve is obtained, but no component name provided to "
          << "calculate reactions from analysis. Use -reaction option";
      exit(1);
    }
    refCurve = readRefCurveData(options::refCurveFilename);
    for (size_t i = 0; i < refCurve.size(); i++) {
      LOG(INFO) << "refCurve[" << i << "] = " << refCurve[i];
    }
  }

  Timer pre_solve(true);
  FEStorage storage;
  NonlinearFESolver solver;
  MeshData md;
  if (!readCdbFile (options::modelFilename, md)) {
    LOG(FATAL) << "Can't read FE info from " << options::modelFilename << "file. exiting..";
  }
  md.compressNumbers();

  Material* mat = CHECK_NOTNULL (MaterialFactory::createMaterial(options::materialName));
  if (mat->getNumC() != options::materialConstants.size())
    LOG(ERROR) << "Material " << options::materialName << " needs exactly " << mat->getNumC()
      << " constants (" << options::materialConstants.size() << " were provided)";
  for (uint16 i = 0; i < mat->getNumC(); i++) {
    mat->Ci(i) = options::materialConstants[i];
  }
  storage.material = mat;
  LOG(INFO) << "Material: " << mat->toString();

  LOG(INFO) << "Loaded components:";
  for (auto& v : md.feComps) {
    LOG(INFO) << v.second;
  }

  // add nodes
  auto sind = storage.createNodes(md.nodesNumbers.size());
  auto ind = md.nodesNumbers;
  for (uint32 i = 0; i < sind.size(); i++) {
    storage.getNode(sind[i]).pos = md.nodesPos[i];
  }


  // add elements
  ind = md.getCellsByAttribute("TYPE", 1);
  sind = storage.createElements(ind.size(), options::elementType);
  for (uint32 i = 0; i < sind.size(); i++) {
    Element& el = storage.getElement(sind[i]);
    assert(el.getNNodes() == md.cellNodes[ind[i]].size());
    for (uint16 j = 0; j < el.getNNodes(); j++) {
      el.getNodeNumber(j) = md.cellNodes[ind[i]][j];
    }
  }

  // add Mpcs
  for (auto& mpc : md.mpcs) {
    // TODO: here the allocated Mpc instance is kept in MeshData,
    // but storage will only have a pointer on it
    storage.addMpc(mpc);
  }

  // add loadBc
  for (auto& v : md.loadBcs) {
    solver.addLoad(v.node, v.node_dof, v.value);
  }

  // add fixBc
  for (auto& v : md.fixBcs) {
    solver.addFix(v.node, v.node_dof, v.value);
  }

  solver.attachFEStorage (&storage);
  solver.numberOfIterations = options::numberOfIterations;
  solver.numberOfLoadsteps = options::numberOfLoadsteps;
    // NOTE: use PARDISO eq. solver by default (if accessible..)
#ifdef NLA3D_USE_MKL
    math::PARDISO_equationSolver eqSolver = math::PARDISO_equationSolver();
    solver.attachEquationSolver(&eqSolver);
#endif


  if (options::useVtk) {
    //obtain job name from path of a FE model file
    std::string jobname = getFileNameFromPath(options::modelFilename);
    VtkProcessor* vtk = new VtkProcessor (&storage, jobname);
    solver.addPostProcessor(vtk);
    vtk->writeAllResults();
  }

  reactProc = NULL;
  if (options::reactionComponentName.length() > 0) {
    reactProc = new ReactionProcessor(&storage);
    solver.addPostProcessor(reactProc);
    FEComponent* feComp = &(md.feComps[options::reactionComponentName]);
    if (feComp->type != FEComponent::NODES) {
      LOG(ERROR) << "To calculate reactions it's needed to provide a component with nodes";
      exit(1);
    }

    LOG(INFO) << "Component for measuring a reaction: " << options::reactionComponentName;
    for (size_t i = 0; i < feComp->list.size(); i++) {
      reactProc->nodes.push_back(feComp->list[i]);
    }
    for (size_t d = 0; d < options::reactionDofs.size(); d++) {
      reactProc->dofs.push_back(options::reactionDofs[d]);
    }
  }

  if (options::rigidBodyMasterNode > 0 && options::rigidBodySlavesComponent.length() > 0) {
    RigidBodyMpc* mpc = new RigidBodyMpc;
    mpc->storage = &storage;
    mpc->masterNode = options::rigidBodyMasterNode;
    FEComponent* comp = &(md.feComps[options::rigidBodySlavesComponent]);
    if (comp->type != FEComponent::NODES) {
      LOG(ERROR) << "For creation a rigid body MPC it's needed to provide a component with nodes";
      exit(1);
    }
    mpc->slaveNodes.assign(comp->list.size(), 0);
    for (size_t i = 0; i < comp->list.size(); i++) {
      mpc->slaveNodes[i] = comp->list[i];
    }

    for (size_t i = 0; i < options::rigidBodyDofs.size(); i++) {
      mpc->dofs.push_back(options::rigidBodyDofs[i]);
    }

    storage.addMpcCollection(mpc);

    std::vector<std::string> slaveDofLabels;
    slaveDofLabels.assign(mpc->dofs.size(), "");
    for (uint16 i = 0; i < mpc->dofs.size(); i++) {
      slaveDofLabels[i] = Dof::dofTypeLabels[i];
    }

    LOG(INFO) << "Rigid Body MPC collection: master node = " << mpc->masterNode << ", number of slaves = " <<  
        mpc->slaveNodes.size() << ", slaves DoFs = " << slaveDofLabels; 
  }

  solver.solve();

  if (reactProc) {
    std::stringstream ss;
    ss << "Reactions:" << std::endl;
    ss << "Loadstep";
    for (size_t d = 0; d < reactProc->dofs.size(); d++) {
      ss << '\t' << Dof::dofTypeLabels[reactProc->dofs[d]];
    }
    ss << std::endl;
    for (size_t ls = 0; ls < options::numberOfLoadsteps+1; ls++) {
      ss << ls;
      for (size_t d = 0; d < reactProc->dofs.size(); d++) {
        ss << '\t' << reactProc->getReactions(reactProc->dofs[d])[ls];
      }
      ss << std::endl;
    }
    LOG(INFO) << ss.str();
    if (refCurve.size() > 0) {
      curCurve = reactProc->getReactions(reactProc->dofs[0]); 
      double _error = compareCurves (refCurve, curCurve);
      LOG(INFO) << "Error between reference loading curve and current is " << _error;
      if (options::curveCompareThreshold > 0.0) {
        if (_error > options::curveCompareThreshold) {
          LOG(ERROR) << "To big error! (upper bound is " << options::curveCompareThreshold << ")";
          exit(1);
        } else {
          LOG(INFO) << "Error is less that the threshold. " << _error << " < " << options::curveCompareThreshold;
        }
      }
    }
  }

  return 0;
}


double compareCurves (const std::vector<double>& refCurve, const std::vector<double>& curCurve) {
  if (refCurve.size() != curCurve.size()) {
    LOG(ERROR) << "compareCurves: tabular data size should be the same for reference curve "
        << "and for current curve. (got refCurve.size() = " << refCurve.size() <<", curCurve.size() = " << curCurve.size() << ")";
    exit(1);
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
    LOG(ERROR) << "readRefCurveData: Expecting third column to a force column";
    exit(1);
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

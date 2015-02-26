#include <conio.h>
#include "sys.h"
#include "FE_Storage.h"
#include "vtk_proc.h"
#include "element_MIXED_8N_3D_P0.h"
#include "Solution.h"
#include "Reaction_proc.h"

void prepeare_processors (FE_Storage_Interface& storage);

// Here is defaults values for command line options
// For command line format see usage() below
uint16 it_num_default = 15;
uint16 ls_num_default = 10;
string mat_name_default = "Compressible Neo-Hookean";
Element::elTypes elTypeDefault = Element::SOLID81;
//TODO: how to initialize vector<double> whithin declaration..
vector<double> mat_Ci_default;
bool is_vtk_default = true;

bool parse_args (int argc, char* argv[], string& model_filename, uint16& it_num, uint16& ls_num, string& mat_name, vector<double>& mat_Ci, bool& is_vtk, Element::elTypes& elType) {
  if (argc < 2) {
    return false;
  }
  model_filename = argv[1];
  argv = argv + 2;
  argc = argc - 2;
  if(cmdOptionExists(argv, argv+argc, "-novtk")) {
    is_vtk = false;
  } else {
    is_vtk = is_vtk_default;
  }

  char* tmp = getCmdOption(argv, argv + argc, "-iterations");
  if (tmp) {
    it_num = atoi(tmp);
  } else {
    it_num = it_num_default;
  }

  tmp = getCmdOption(argv, argv + argc, "-loadsteps");
  if (tmp) {
    ls_num = atoi(tmp);
  } else {
    ls_num = ls_num_default;
  }

  tmp = getCmdOption(argv, argv + argc, "-element");
  if (tmp) {
    elType = Element::elName2elType(tmp);
  } else {
    elType = elTypeDefault;
  }

  vector<char*> vtmp = getCmdManyOptions(argv, argv + argc, "-material");
  if (vtmp.size() == 0) {
    mat_name = mat_name_default;
    //TODO: how to initialize vector<double> whithin declaration..
    mat_Ci.push_back(30.0);
    mat_Ci.push_back(0.499);
    //mat_Ci = mat_Ci_default;
  } else {
    mat_name = vtmp[0];
    for (uint16 i = 1; i < vtmp.size(); i++) {
      mat_Ci.push_back(atof(vtmp[i]));
    }
  }

  return true;
}

void usage () {
  echolog("nla3d model_file [-material mat_name mat_C0 mat_C1 ..] [-iterations it_num] [-loadsteps ls_num] [-novtk]");
}

int main (int argc, char* argv[])
{
  uint16 it_num;
  uint16 ls_num;
  string mat_name;
  string model_filename;
  bool is_vtk;
  vector<double> mat_Ci;
  Element::elTypes elType;

	echolog("---=== WELCOME TO NLA PROGRAM ===---");
  if (!parse_args(argc, argv, model_filename, it_num, ls_num, mat_name, mat_Ci, is_vtk, elType)) {
    usage();
    exit(1);
  }
	Timer pre_solve(true);
	FE_Storage storage;
  storage.elType = elType;
  if (!read_ans_data(model_filename.c_str(), &storage)) {
    error("Can't read FE info from %s file. exiting..", model_filename.c_str());
  }
  Material* mat = createMaterial(mat_name);
  if (mat->getNumC() != mat_Ci.size())
    error("Material %s needs exactly %d constants (%d were provided)", mat_name.c_str(), mat->getNumC(), mat_Ci.size());
  for (uint16 i = 0; i < mat->getNumC(); i++) 
    mat->Ci(i) = mat_Ci[i];
	storage.material = mat;

	Solution sol;
	sol.attach(&storage);
	sol.setqIterat(it_num);
	sol.setqLoadstep(ls_num);
  if (is_vtk) {
    char _jobname[100];
    char dummy[1000];
    //obtain job name from path of a FE model file
    //_splitpath is OS dependent (MS VS)
    _splitpath(model_filename.c_str(), dummy, dummy, _jobname, dummy);
    string jobname(_jobname);
    Vtk_proc* vtk = new Vtk_proc(&storage, jobname);
  }
  prepeare_processors(storage); 

	echolog("Preprocessor time: %f sec.", pre_solve.stop());
	sol.run();
	return 0;
}

//while performing an analysis we need to process and store some results data
//here is a mechanism of processors to deal with it
//but it needs a lot of options to setup it
//now it's needed to hardcode initialization and setup of a such processors
void prepeare_processors (FE_Storage_Interface& storage) {
	Reaction_proc* proc = new Reaction_proc(&storage, "loading_force.txt");
	list<BC_dof_constraint> &lst = storage.get_dof_const_BC_list();

	list<BC_dof_constraint>::iterator bc_dof = lst.begin();
	while (bc_dof != lst.end())
	{
		if (fabs(bc_dof->value) > 0)
		{
			proc->nodes.push_back(bc_dof->node);
			proc->dofs.push_back(bc_dof->node_dof);
		}
		bc_dof++;
	}
}


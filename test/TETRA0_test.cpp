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
#include "elements/TETRA0.h"
#include <tuple>

using namespace nla3d;

typedef std::vector<std::tuple<double, double, double>> disp_vec_t;
typedef std::vector<math::MatSym<3>> stress_vec_t;

disp_vec_t readDispData (std::string filename);
stress_vec_t readStressData (std::string filename);

int main (int argc, char* argv[]) {
    std::string cdb_filename;
    std::string res_disp_filename;
    std::string res_stress_filename;

    if (argc > 1) {
      cdb_filename = argv[1];
    } else {
      LOG(FATAL) << "You should provide the path to mesh (cdb file)";
    }
    if (argc > 2) {
      res_disp_filename = argv[2];
    }
    if (argc > 3) {
      res_stress_filename = argv[3];
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
    sind = storage.createElements(ind.size(), ElementType::TETRA0);
    for (uint32 i = 0; i < sind.size(); i++) {
      ElementTETRA0& el = dynamic_cast<ElementTETRA0&>(storage.getElement(sind[i]));
      el.getNodeNumber(0) = md.cellNodes[ind[i]][0];
      el.getNodeNumber(1) = md.cellNodes[ind[i]][1];
      el.getNodeNumber(2) = md.cellNodes[ind[i]][2];
      el.getNodeNumber(3) = md.cellNodes[ind[i]][4];
      el.E = 1.0e8;
      el.my = 0.3;
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
	solver.attachFEStorage (&storage);

	// VtkProcessor* vtk = new VtkProcessor (&storage, "tetra");
    // solver.addPostProcessor(vtk);
    // vtk->writeAllResults();

    solver.solve();

    // check displacements results with Ansys data
    if (res_disp_filename != "") {
        auto ans_disps = readDispData(res_disp_filename);
        for (uint32 i = 1; i <= storage.nNodes(); i++) {
            const double dx = storage.getNodeDofSolution(i, Dof::UX);
            const double dy = storage.getNodeDofSolution(i, Dof::UY);
            const double dz = storage.getNodeDofSolution(i, Dof::UZ);
            const auto& ans_disp = ans_disps[i - 1];
            CHECK_EQTH(dx, std::get<0>(ans_disp), 1.0e-10);
            CHECK_EQTH(dy, std::get<1>(ans_disp), 1.0e-10);
            CHECK_EQTH(dz, std::get<2>(ans_disp), 1.0e-10);
        }
    }

    // check stress results with Ansys data
    if (res_stress_filename != "") {
        auto ans_stresses = readStressData(res_stress_filename);
        math::MatSym<3> mat;
        for (uint32 i = 1; i <= storage.nNodes(); i++) {
            mat.zero();
            storage.getElement(i).getTensor(&mat, tensorQuery::E);
            CHECK(mat.compare(ans_stresses[i - 1], 1.0e-3));
        }
    }
}

disp_vec_t readDispData (std::string filename) {
    disp_vec_t disp_vec;
    std::ifstream file(filename);
    if (!file.is_open()) {
        return disp_vec;
    }
    double tmp, dx, dy, dz;
    // file format:
    //  node UX   UY   UY
    //        1.  0.00000000E+00  0.00000000E+00  0.00000000E+00
    //        2. -0.34498949E-06 -0.32502630E-06  0.00000000E+00
    //  ...

    // skip first line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (!file.eof()) {
      file >> tmp >> dx >> dy >> dz;
      disp_vec.emplace_back(dx, dy, dz);
    }
    return disp_vec;
}

stress_vec_t readStressData (std::string filename) {
    stress_vec_t stress_vec;
    std::ifstream file(filename);
    if (!file.is_open()) {
        return stress_vec;
    }
    double tmp, sx, sy, sz, sxy, syz, sxz;
    // file format:
    //    elem SX   SY   SZ   SXY  SYZ  SXZ
    //          1.  0.28688422E+02  0.90668364E+01  0.11326577E+02  0.65410609E-01 -0.36082213E+01 -0.19945953E+01
    //          2.  0.24542393E+02  0.10518168E+02  0.10518168E+02 -0.23288821E+01  0.00000000E+00  0.00000000E+00
    //  ...

    // skip first line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    math::MatSym<3> mat;
    while (!file.eof()) {
      file >> tmp >> mat.comp(0, 0) >> mat.comp(1, 1) >> mat.comp(2, 2)
                  >> mat.comp(0, 1) >> mat.comp(1, 2) >> mat.comp(0, 2);
      stress_vec.push_back(mat);
    }
    return stress_vec;
}

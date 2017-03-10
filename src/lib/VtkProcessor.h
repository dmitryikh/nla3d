// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include <set>
#include <map>
#include "sys.h"
#include "PostProcessor.h"
#include "FEStorage.h"
#include "elements/query.h"


namespace nla3d {

// VtkProcessor is dedicated to dump mesh and results into *.vtk legacy textual format. On pre()
// callback filename0.vtk is written which contains initial state of FE problem before solution. On
// every solution step process(...) callback is called which writes filename%d.vtk with solution
// data relevant to this step.  Later this banch of files could be used to visualize solution fields
// in external program (ex.  Paraview).
//
// By default, VtkProcessor write DoF fields (displacements, reactions for structural problems and
// temperatures for thermal). But by registering element results by registerResult(...) one can add
// addition fields (scalar, vector or tensor ones) to be written into vtk files.

class VtkProcessor : public PostProcessor {
public:
	VtkProcessor(FEStorage *st, std::string _fileName);
	virtual ~VtkProcessor();
	virtual void pre();
	virtual void process(uint16 curLoadstep);
	virtual void post(uint16 curLoadstep);

  // One can use writeAllResults(true) in order to ask VtkProcess to determine which query codes is
  // relevant to FEStorage's elements automatically. As results, all accessible element results will
  // be written into vtk files.
  void writeAllResults(bool write = true);

  // The methods below can be used to manually set query codes for elements. Every query will be add
  // addition cell data fields into vtk file.
  bool registerResults(std::initializer_list<scalarQuery> queries);
  bool registerResults(std::initializer_list<vectorQuery> queries);
  bool registerResults(std::initializer_list<tensorQuery> queries);

private:
	std::string file_name;

	void write_header(std::ofstream &file);
  // write geometry(mesh) into the vtk's `file`. If `def == true` then write node coordinates in
  // deformed state
	void write_geometry(std::ofstream &file, bool def = false);
	void write_point_data(std::ofstream &file);
	void write_cell_data(std::ofstream &file);

  void writeScalar(std::ofstream &file, const char* name, std::vector<double>& data);
  void writeTensor(std::ofstream &file, const char* name, std::vector<math::MatSym<3> >& data);
  void writeVector(std::ofstream &file, const char* name, std::vector<math::Vec<3> >& data);

  void revealAllResults();

  bool useDisplacementsDofs = false;
  bool _writeAllResults = false;

  std::set<Dof::dofType> nodalDofTypes;
  std::set<Dof::dofType> elementDofTypes;

  std::set<scalarQuery> cellScalarQueries;
  std::set<vectorQuery> cellVectorQueries;
  std::set<tensorQuery> cellTensorQueries;
};

} // namespace nla3d

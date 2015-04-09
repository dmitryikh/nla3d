// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "PostProcessor.h"
#include "FEStorage.h"
#include "elements/query.h"

namespace nla3d {

class VtkProcessor : public PostProcessor {
public:
	VtkProcessor(FEStorage *st, std::string _fileName);
	virtual ~VtkProcessor();
	virtual void pre ();
	virtual void process (uint16 curLoadstep);
	virtual void post (uint16 curLoadstep);
private:
	std::string file_name;

	void write_header (std::ofstream &file);
	void write_geometry(std::ofstream &file, bool def = false);
	void write_point_data(std::ofstream &file);
	void write_cell_data(std::ofstream &file);

  void writeScalar(std::ofstream &file, const char* name, std::vector<double>& data);
  void writeTensor(std::ofstream &file, const char* name, std::vector<math::MatSym<3> >& data);
  void writeVector(std::ofstream &file, const char* name, std::vector<math::Vec<3> >& data);

  std::vector<query::scalarQuery> cellScalarQuery;
  std::vector<query::vectorQuery> cellVectorQuery;
  std::vector<query::tensorQuery> cellTensorQuery;
private:
};

} // namespace nla3d

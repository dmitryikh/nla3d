// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "PostProcessor.h"
#include "FEStorage.h"

namespace nla3d {

class ReactionProcessor : public PostProcessor {
public:
	ReactionProcessor(FEStorage *st);
	ReactionProcessor(FEStorage *st, std::string _filename);
	virtual ~ReactionProcessor() { };

	virtual void pre ();
	virtual void process (uint16 curLoadstep);
	virtual void post (uint16 curLoadstep);

  std::vector<double> getReactions (Dof::dofType dof);
	
	std::vector<uint32> nodes;
	std::vector<Dof::dofType> dofs;
protected:
  std::string filename;
	std::vector<std::vector<double> > sumOfDofsReactions;
};

} // namespace nla3d

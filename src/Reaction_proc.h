// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "PostProcessor.h"
#include "FEStorage.h"

namespace nla3d {

class Reaction_proc : public PostProcessor {
public:
	Reaction_proc(FEStorage *st);
	Reaction_proc(FEStorage *st, std::string _filename);
	virtual ~Reaction_proc() { };
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
  std::vector<double> getReactions ();
	
	std::vector<uint32> nodes;
	std::vector<uint16> dofs;
protected:
  std::string filename;
	std::vector<double> reactVec;
};

} // namespace nla3d

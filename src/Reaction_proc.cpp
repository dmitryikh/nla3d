// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "Reaction_proc.h"

namespace nla3d {

Reaction_proc::Reaction_proc(FEStorage *st) : PostProcessor(st)
{
	name ="Reaction_proc";
	active = true;
  reactVec.clear();
}

Reaction_proc::Reaction_proc(FEStorage *st, std::string _filename) : PostProcessor(st)
{
	name ="Reaction_proc";
	active = true;
  reactVec.clear();
  filename = _filename;
}

void Reaction_proc::pre(uint16 qLoadstep)
{
	if (nodes.size() == 0 || dofs.size() == 0 || dofs.size() != nodes.size()) {
    warning("Reaction_proc::pre: can't work. No nodes and/or dofs. Processor name = %s", name.c_str()); 
		active = false;
  }
	if (!active) return;
  if (filename.length() > 0) {
    std::ofstream file(filename.c_str(),std::ios::trunc);
    if (!file)
    {
      warning("Reaction_proc::pre: Can't create a file with name %s", filename.c_str());
      return;
    }
    file << 0.0 << std::endl;
    file.close();
  }

  reactVec.resize(qLoadstep+1);
  reactVec[0] = 0.0;
}

void Reaction_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	double force = 0.0;
	for (uint32 i=0; i < nodes.size(); i++)
		force += storage->get_reaction_force(nodes[i],dofs[i]);
  if (filename.length() > 0) {
    std::ofstream file(filename.c_str(),std::ios::app);
    if (!file)
    {
      warning("Reaction_proc::pre: Can't create a file with name %s", filename.c_str());
      return;
    }
    file << force << std::endl;
    file.close();
  }
  reactVec[curLoadstep] = force;
}


void Reaction_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
}


std::vector<double>  Reaction_proc::getReactions () {
  return reactVec;
}

} // namespace nla3d

#include "Reaction_proc.h"

Reaction_proc::Reaction_proc(FE_Storage_Interface *st) : Post_proc(st)
{
	name ="Reaction_proc";
	active = true;
  reactVec.clear();
}

Reaction_proc::Reaction_proc(FE_Storage_Interface *st, string _filename) : Post_proc(st)
{
	name ="Reaction_proc";
	active = true;
  reactVec.clear();
  filename = _filename;
}

void Reaction_proc::pre(uint16 qLoadstep)
{
	if (nodes.size() == 0 || dofs.size() == 0 || dofs.size() != nodes.size()) {
    warning("Reaction_proc::pre: can't work. No nodes and/or dofs. Processor name = %s", name); 
		active = false;
  }
	if (!active) return;
  if (filename.length() > 0) {
    ofstream file(filename,ios::trunc);
    if (!file)
    {
      warning("Reaction_proc::pre: Can't create a file with name %s", filename);
      return;
    }
    file << 0.0 << endl;
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
    ofstream file(filename,ios::app);
    if (!file)
    {
      warning("Reaction_proc::pre: Can't create a file with name %s", filename);
      return;
    }
    file << force << endl;
    file.close();
  }
  reactVec[curLoadstep] = force;
}


void Reaction_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
}


vector<double>  Reaction_proc::getReactions () {
  return reactVec;
}

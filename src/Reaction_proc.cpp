#include "Reaction_proc.h"

Reaction_proc::Reaction_proc(FE_Storage_Interface *st) : Post_proc(st)
{
	name ="Reaction_proc";
	active = true;
}

void Reaction_proc::pre(uint16 qLoadstep)
{
	if (nodes.size() == 0 && dofs.size() == 0 && dofs.size() != nodes.size())
		active = false;
	if (!active) return;
	stringstream ss;
	ss << "POST" << nPost_proc << "_" << name << ".txt";
	filename = ss.str();
	ofstream file(filename,ios::trunc);
	if (!file)
	{
		warning("Reaction_proc::pre: Can't create a file with name %s", filename);
		return;
	}
	file << 0.0 << endl;
	file.close();
}

void Reaction_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	double force = 0.0;
	for (uint32 i=0; i < nodes.size(); i++)
		force += storage->get_reaction_force(nodes[i],dofs[i]);
	ofstream file(filename,ios::app);
	if (!file)
	{
		warning("Reaction_proc::process: Can't create a file with name %s", filename);
		return;
	}
	file << force << endl;
	file.close();
}


void Reaction_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
}

//void Reaction_proc::load_from_file (const char *filename)
//{
//	ifstream file(filename);
//	if (!file)
//	{
//		warning("RForce_proc::load_from_file: Can't open input file `%s`",filename);
//		return;
//	}
//	char buf[1024]="";
//	uint32 val;
//	while (file.getline(buf, 1024))
//	{
//		sscanf(buf, "%d", &val);
//		if (val < 1) break;
//		nodes.push_back(val);
//	}
//	file.close();
//}
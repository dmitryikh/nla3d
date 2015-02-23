#include "Cmd_Shell.h"

Cmd_Shell::Cmd_Shell(FE_Storage_Interface *_storage)
{
	solution.attach(storage);
	storage = _storage;
	Command cmd;
	cmd.describtion="Provide the user with the list of the commands and outline short description for each.";
	cmd.shortdescribption="List the NLA's commands";
	cmd.name="help";
	cmd.p_func=help_command;
	commands.push_back(cmd);
	cmd.describtion="Show large information about FEM and Solver state (nNodes, nElements, nBCs, nBand, criteria, nLoadStep, analysis type, nInt, Material type, Material properties, maxIterat).";
	cmd.shortdescribption="Show status of FEM and Solver";
	cmd.name ="status";
	cmd.p_func=status_command;
	commands.push_back(cmd);
	//cmd.describtion="Load FEM data from file. Old data would be replaced. Usage: loadmodel(filename.mdl)";
	//cmd.shortdescribption="Load FEM data from mdl file";
	//cmd.name="loadmodel";
	//cmd.p_func =loadmodel_command;
	//commands.push_back(cmd);
	//cmd.describtion="Save FEM data to file. Solution options isn't written. Usage: savemodel(filename.mdl)";
	//cmd.shortdescribption="Save FEM data to mdl file";
	//cmd.name="savemodel";
	//cmd.p_func =savemodel_command;
	//commands.push_back(cmd);
	//cmd.describtion="Set material for analysis. Usage: material(material_type, C1,C2)\nCompressible Neo-Hookean material: material(Comp_Neo_Hookean, E, mu)\nSimple Hookean model material: material(Hookean, E, mu)\n";
	//cmd.shortdescribption="Set material for analysis";
	//cmd.name="mat";
	//cmd.p_func =mat_command;
	//commands.push_back(cmd);
	cmd.describtion="Set convergence criteria for solver. Usage: criteria(1.0e-5)";
	cmd.shortdescribption="Set convergence criteria for solver";
	cmd.name="crit";
	cmd.p_func =crit_command;
	commands.push_back(cmd);
	//cmd.describtion="Set analysis type. Usage: analysistype(type)\ntype: AXYM - axymmetric, PLANE - plane strain";
	//cmd.shortdescribption="Set analysis type";
	//cmd.name="antype";
	//cmd.p_func = antype_command;
	//commands.push_back(cmd);
	cmd.describtion="Set the number of loadsteps for solver. Usage: loadsteps(15)";
	cmd.shortdescribption="Set number of loadsteps for solver";
	cmd.name="steps";
	cmd.p_func = steps_command;
	commands.push_back(cmd);
	cmd.describtion="Set the number of maximum iterations per step for solver. Usage: maxiter(20)";
	cmd.shortdescribption="Set the number of maximum iterations per step for solver";
	cmd.name="iters";
	cmd.p_func =iters_command;
	commands.push_back(cmd);
	cmd.describtion="Set the integration scheme for elements. Usage: intscheme(3)\nSchemes: 1x1\t2x2\t3x3\t4x4";
	cmd.shortdescribption="Set the integration scheme for elements";
	cmd.name="int";
	cmd.p_func =int_command;
	commands.push_back(cmd);
	cmd.describtion="Run the solution process.";
	cmd.shortdescribption="Run the solution process";
	cmd.name="solve";
	cmd.p_func =solve_command;
	commands.push_back(cmd);
	cmd.describtion="Multiply BC's values to x. Usage: multbc(x)";
	cmd.shortdescribption="Multiply BC's values to x";
	cmd.name="multbc";
	cmd.p_func =multbc_command;
	commands.push_back(cmd);
	//cmd.describtion="Make strip model. Usage: maketestmodel(B,H,Bn,Hn,alpha,r,delta)";
	//cmd.shortdescribption="Make strip model";
	//cmd.name="stripmodel";
	//cmd.p_func =stripmodel_command;
	//commands.push_back(cmd);

	//cmd.describtion="read command input file";
	//cmd.shortdescribption="read command input file";
	//cmd.name="input";
	//cmd.p_func = input_command;
	//commands.push_back(cmd);

	cmd.describtion="Load a FEM model from the ansys *.cdb file. Usage: ansload(filename)";
	cmd.shortdescribption="Load a FEM model from the ansys *.cdb file";
	cmd.name="ansload";
	cmd.p_func = ansload_command;
	commands.push_back(cmd);

	//cmd.describtion="adding egraph_proc to solution process. Usage: egraph_proc(EX). valid input: EX,EY,EZ,EXY,EXZ,EYZ,EVOL,SX,SY,SZ and so on";
	//cmd.shortdescribption="adding egraph_proc to solution process";
	//cmd.name="egraph_proc";
	//cmd.p_func = egraph_proc_command;
	//commands.push_back(cmd);

	//cmd.describtion="Make circ model. Usage: circmodel(r1,r2,fi1,fi2,nR,nFi,P)";
	//cmd.shortdescribption="Make circ model";
	//cmd.name="circmodel";
	//cmd.p_func = circmodel_command;
	//commands.push_back(cmd);
}

void Cmd_Shell::run()
{
	greeting();
	command_processor();
	//bye();
}
void Cmd_Shell::command_processor()
{
	string command;
	char buff[512];
	while (command != "exit")
	{
		cout << '>';
		cin.getline(buff,512);
		log(">%s",buff);
		command.assign(buff);
		atomize(command);
		if (atoms.size() == 0) continue;
		bool find=false;
		for (uint16 i=0; i<commands.size() && find==false;i++)
		{
			if (commands[i].name == atoms[0])
			{
				commands[i].p_func(*this);
				find=true;
			}
		}
		if (!find && atoms[0] != "exit")
			warning("Unknown command. Please, re-type");
	}
}
void Cmd_Shell::atomize (string cmd)
{
	char delimeters[] = "(),";
	atoms.clear();
	if (cmd.length()==0) return;
	uint16 start=0;
	uint16 i = 0;
	for (i=0; i < cmd.length(); i++)
	{
		for (char *p=delimeters; *p != 0; p++)
		{
			if (cmd[i]==*p)
			{
				atoms.push_back(string(cmd,start,i-start));
				atoms.push_back(string(cmd,i,1));
				start = i+1;
			}
		}
	}
	if (start!=i)
		atoms.push_back(string(cmd,start,i-start));
	//удалим пробелы
	for (uint16 j=0;j<atoms.size();j++)
		del_spaces(atoms[j]);
}

void help_command(Cmd_Shell &shell)
{
	if (shell.atoms.size()==1)
	{
		for (uint16 i=0;i<shell.commands.size(); i++)
			echolog("%s - %s", shell.commands[i].name.c_str(), shell.commands[i].shortdescribption.c_str());
		return;
	}
	else if (shell.atoms[1]=="(" && shell.atoms.size()>2)
	{
		bool find = false;
		for (uint16 i=0;i<shell.commands.size();i++)
		{
			if (shell.commands[i].name == shell.atoms[2])
			{
				echolog("%s - %s",shell.commands[i].name.c_str(), shell.commands[i].describtion.c_str());
				find=true;
			}
		}
		if (find==false)
		{
			warning("command %s isn't found.", shell.atoms[2].c_str());
			return;
		}
		return;
	}
	warning ("mistakes in command.");
}

void status_command (Cmd_Shell &shell)
{
	if (shell.atoms.size()>1)
	{
		warning("mistakes in command.");
		return;
	}
	echolog("%d nodes,\t%d elements,\t%d B.C.'s",shell.storage->getNumNode(), shell.storage->getNumElement(), shell.storage->getNumBound());
	echolog("nBand = %d,\tqLoadStep = %d", shell.storage->get_nBand(), shell.solution.getqLoadstep());
	echolog("maxIterat = %d,\tcriteria = %f",shell.solution.getqIterat(), shell.solution.getCriteria());
	//TODO: make good class for material
	//if (shell.storage->getMaterial() != NULL)
	//	echolog("Material: %s\tE = %f,\tmu = %f",shell.storage->getMaterial().getName(), shell.storage.material->E, shell.storage.material->mu);
	//else
	//	echolog("Material: not defined");
}

//void loadmodel_command (Cmd_Shell &shell)
//{
//	if (shell.atoms.size() != 4)
//	{
//		warning("mistakes in command.");
//		return;
//	}
//	if(shell.storage.load_mdl(shell.atoms[2].c_str()))
//		echolog("%s loaded: %d nodes, %d elements, %d bcs", shell.atoms[2].c_str(), shell.storage.nn, shell.storage.en, shell.storage.bounds.size());
//	else
//		warning("%s wasn't loaded due to errors", shell.atoms[2].c_str());
//}

//void savemodel_command (Cmd_Shell &shell)
//{
//	if (shell.atoms.size() != 4)
//	{
//		warning("mistakes in command.");
//		return;
//	}
//	if(shell.storage.save_mdl(shell.atoms[2].c_str()))
//		echolog("FEM data saved", shell.atoms[2].c_str(), shell.storage.nn, shell.storage.en, shell.storage.bounds.size());
//	else
//		warning("Data wasn't saved to file '%s' due to errors", shell.atoms[2].c_str());
//}


//TODO: make good mat class
//void mat_command (Cmd_Shell &shell)
//{
//	if (shell.atoms.size() != 8)
//	{
//		warning("mistakes in command.");
//		return;
//	}
//	if (shell.storage.material)
//		delete shell.storage.material;
//	if (shell.atoms[2] == "Comp_Neo_Hookean")
//		shell.storage.material = new Material_Comp_Neo_Hookean();
//	else if (shell.atoms[2] == "Hookean")
//		shell.storage.material = new Material_Hookean();
//	else
//	{
//		warning("Unknown type of material");
//		return;
//	}
//	shell.storage.material->E = atof(shell.atoms[4].c_str());
//	shell.storage.material->mu = atof(shell.atoms[6].c_str());
//}

void crit_command (Cmd_Shell &shell)
{
	if (shell.atoms.size() != 4)
	{
		warning("mistakes in command.");
		return;
	}
	shell.solution.setCriteria(atof(shell.atoms[2].c_str()));
}

//void antype_command (Cmd_Shell &shell)
//{
//	if (shell.atoms.size() != 4)
//	{
//		warning("mistakes in command.");
//		return;
//	}
//	if (shell.atoms[2] == "MIXED")
//		shell.storage.fem_type = FEM_MIXED;
//	/*else if (shell.atoms[2] == "AXYM")
//		shell.storage.fem_type = FEM_AXYM;
//	else if (shell.atoms[2] == "PLANE_MIXED")
//		shell.storage.fem_type = FEM_PLANE_MIXED;*/
//	else
//	{
//		warning("Unknown type of analysis");
//		return;
//	}
//}

void steps_command (Cmd_Shell &shell)
{
	if (shell.atoms.size() != 4)
	{
		warning("mistakes in command.");
		return;
	}
	shell.solution.setqLoadstep(atoi(shell.atoms[2].c_str()));
}

void iters_command (Cmd_Shell &shell)
{
	if (shell.atoms.size() != 4)
	{
		warning("mistakes in command.");
		return;
	}
	shell.solution.setqIterat(atoi(shell.atoms[2].c_str()));
}


void int_command (Cmd_Shell &shell)
{
	if (shell.atoms.size() != 4)
	{
		warning("mistakes in command.");
		return;
	}
	uint16 sch = atoi(shell.atoms[2].c_str());
	if (sch < 1 || sch > 4)
	{
		warning("wrong integration scheme.");
		return;
	}
	Element::set_nInt(sch);
}

void solve_command (Cmd_Shell &shell)
{
	if (shell.atoms.size() != 1)
	{
		warning("mistakes in command.");
		return;
	}
	shell.solution.run();
}

void multbc_command (Cmd_Shell &shell)
{
	if (shell.atoms.size() != 4)
	{
		warning("mistakes in command.");
		return;
	}
	double mult = atof(shell.atoms[2].c_str());
	for (uint32 i=0; i<shell.storage->getNumBound(); i++)
		shell.storage->getBound(i).value *= mult;
}

//void stripmodel_command(Cmd_Shell &shell)
//{
//	if (shell.atoms.size() != 16)
//	{
//		warning("mistakes in command.");
//		return;
//	}
//	shell.storage.makestripmodel(atof(shell.atoms[2].c_str()), atof(shell.atoms[4].c_str()),
//		atoi(shell.atoms[6].c_str()), atoi(shell.atoms[8].c_str()), atof(shell.atoms[10].c_str())*M_PI/180.0f,atof(shell.atoms[12].c_str()),atof(shell.atoms[14].c_str()));
//	Def_Form_proc *def_proc = new Def_Form_proc(&shell.storage);
//}

void input_command (Cmd_Shell &shell)
{
	if (shell.atoms.size() != 4)
	{
		warning("mistakes in command.");
		return;
	}
	char buf[1024];
	ifstream file;
	string command;
	file.open(shell.atoms[2].c_str());
	if (!file)
	{
		warning("Can't open file %s", shell.atoms[2].c_str());
		return;
	}
	while (file.getline(buf,1024))
	{
		command.assign(buf);
		shell.atomize(command);
		if (shell.atoms.size() == 0) continue;
		bool find=false;
		for (uint16 i=0; i<shell.commands.size() && find==false;i++)
		{
			if (shell.commands[i].name == shell.atoms[0])
			{
				shell.commands[i].p_func(shell);
				find=true;
			}
		}
		if (!find && shell.atoms[0] != "exit")
			warning("Unknown command. Please, re-type");
	}
}

void ansload_command (Cmd_Shell &shell)
{
	if (shell.atoms.size() != 4)
	{
		warning("mistakes in command.");
		return;
	}
	if(read_ans_data(shell.atoms[2].c_str(), shell.storage))
	{
		echolog("%s loaded: %d nodes, %d elements, %d bcs", shell.atoms[2].c_str(), shell.storage->getNumNode(), shell.storage->getNumElement(), shell.storage->getNumBound());
	}
	else
		warning("%s wasn't loaded due to errors", shell.atoms[2].c_str());
}

//void egraph_proc_command (Cmd_Shell &shell)
//{
//	if (shell.atoms.size() != 4)
//	{
//		warning("mistakes in command.");
//		return;
//	}
//	//TODO: добавить выбор точек интегрирования
//	EGraph_proc * egraph_proc = new EGraph_proc(&shell.storage, shell.atoms[2],4);
//	echolog("EGraph_proc No %d was added", egraph_proc->getProc_num());
//}

//void circmodel_command(Cmd_Shell &shell)
//{
//	if (shell.atoms.size() != 16)
//	{
//		warning("mistakes in command.");
//		return;
//	}
//	shell.storage.makecircmodel(atof(shell.atoms[2].c_str()), atof(shell.atoms[4].c_str()),
//		atof(shell.atoms[6].c_str())*M_PI/180.0f, atof(shell.atoms[8].c_str())*M_PI/180.0f, atoi(shell.atoms[10].c_str()),atoi(shell.atoms[12].c_str()),atof(shell.atoms[14].c_str()));
//	Def_Form_proc *def_proc = new Def_Form_proc(&shell.storage);
//}
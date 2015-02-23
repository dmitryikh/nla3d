#pragma once
#include "FE_Storage.h"
#include "Solution.h"

class Cmd_Shell;

struct Command
{
	string describtion;
	string shortdescribption;
	string name;
	void (*p_func)(Cmd_Shell &shell);
};

void help_command (Cmd_Shell &shell);
void status_command (Cmd_Shell &shell);
//void loadmodel_command (Cmd_Shell &shell);
//void savemodel_command (Cmd_Shell &shell);
//void mat_command (Cmd_Shell &shell);
void crit_command (Cmd_Shell &shell);
//void antype_command (Cmd_Shell &shell);
void steps_command (Cmd_Shell &shell);
void iters_command (Cmd_Shell &shell);
void int_command (Cmd_Shell &shell);
void solve_command (Cmd_Shell &shell);
void multbc_command (Cmd_Shell &shell);
//void stripmodel_command(Cmd_Shell &shell);
//void input_command (Cmd_Shell &shell);
//void egraph_proc_command (Cmd_Shell &shell);
void ansload_command (Cmd_Shell &shell);
//void circmodel_command (Cmd_Shell &shell);

class Cmd_Shell
{
public:
	Cmd_Shell (FE_Storage_Interface *_storage);
	void run ();
	void greeting ()
	{
		echolog("---Welcome to NLA 3D program!---\nVersion: %s\t Data: %s", SYS_VERSION, SYS_DATA);
		echolog("Try 'help (command)' to see command description or 'help' to see all commands");
	}
	void command_processor ();
	void atomize (string cmd);
	FE_Storage_Interface* storage;
	Solution solution;
	vector<string> atoms;
	vector<Command> commands;
	//commands
};

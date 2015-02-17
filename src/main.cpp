#include <conio.h>
#include "sys.h"
#include "math\Mat.h"
#include "math\Mat_Band.h"
#include "FE_Storage.h"
#include "vtk_proc.h"
#include "basic_elements.h"
//#include "Render.h"
#include "Solution.h"
#include <mkl.h>
#include "math\Sparse_Matrix.h"
#include <limits>
#include "Reaction_proc.h"

void main ()
{
	echolog("---=== WELCOME TO NLA PROGRAM ===---");
	Timer pre_solve(true);
	FE_Storage<MIXED_8N_3D_P0_NL> storage;
	read_ans_data("model.cdb", &storage);
	Material_Comp_Neo_Hookean* mat = new Material_Comp_Neo_Hookean();
	//Material_Hookean* mat = new Material_Hookean();
	mat->Ci(0) = 30;
	mat->Ci(1) = 0.499;
	storage.material = mat;
	Vtk_proc* vtk = new Vtk_proc(&storage);
	Solution sol;
	sol.attach(&storage);
	sol.setqIterat(15);
	sol.setqLoadstep(10);
	double a1 = 130*M_PI/180;
	double a2 = 50*M_PI/180;
	ifstream in("LIST.txt");
	uint32 node;
	double cx,cz;
	double fi_sec = 2*M_PI/3;
	cx = cos(fi_sec+M_PI/2);
	cz = -sin(fi_sec+M_PI/2);
	BC_MPC mpc;
	mpc.b = 0.0;
	while (in)
	{	
		mpc.eq.clear();
		node = 0;
		in >> node;
		if (node)
		{
			mpc.eq.push_back(MPC_token(node,0,cx));
			mpc.eq.push_back(MPC_token(node,2,cz));
			//storage.add_bounds(mpc);
		}
	}
	in.close();
	Reaction_proc proc(&storage);
	list<BC_dof_constraint> &lst = storage.get_dof_const_BC_list();
	list<BC_dof_constraint>::iterator bc_dof = lst.begin();
	while (bc_dof != lst.end())
	{
		if (fabs(bc_dof->value) > 0)
		{
			proc.nodes.push_back(bc_dof->node);
			proc.dofs.push_back(bc_dof->node_dof);
		}
		bc_dof++;
	}
	/*BC_MPC mpc1;
	mpc1.b = 0.0;
	mpc1.eq.push_back(MPC_token(5,0,cos(a1)));
	mpc1.eq.push_back(MPC_token(5,1,sin(a1)));
	BC_MPC mpc2;
	mpc2.b = 0.3;
	mpc2.eq.push_back(MPC_token(3,0,cos(a2)));
	mpc2.eq.push_back(MPC_token(3,1,sin(a2)));*/
	//storage.add_bounds(mpc1);
	//storage.add_bounds(mpc2);
	echolog("Preprocessor time: %f sec.", pre_solve.stop());
	sol.run();
	//Mat_test();
	/*Cmd_Shell *shell = new Cmd_Shell();
	shell->run();
	delete shell;
	long t2 = clock();
	echolog("Program is closed (Worktime= %10.2fs)",(float) (t2-t1)/1000.0f);*/
	return;
}




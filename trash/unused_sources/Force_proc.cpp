#include "Force_proc.h"
#include "math\Vec.h"
#include "math\Mat.h"

Force_by_stress_proc::Force_by_stress_proc(FE_Storage *st) : Post_proc(st)
{
	name = "Force_by_stress";
}

void Force_by_stress_proc::pre(uint16 qLoadstep)
{
	/*echolog("STRESS_PROC INFORMATION:\nEl\tn1\tn2");
	for (uint16 i=0; i < elems.size(); i++)
	{
		echolog("%d.\t%d\t%d\t%d",i+1,elems[i],nodes1[i],nodes2[i]);
	}*/
	if (!active) return;
	forces.clear();
	forces.push_back(0.0f);
}

void Force_by_stress_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	// пока реализуем только для ПНС
	if (storage->fem_type != FEM_PLANE)
	{
		warning("Stress_proc::process: this post_proc works only with PLANE analysis");
		return;
	}
	if (_4N_2D_Form::nInt != 3)
	{
		warning("Stress_proc::process: this post_proc works only with integration scheme 3x3");
		return;
	}
	Mat<3,3> matX;
	Mat<3,3> matS;
	Mat<3,3> matT;
	double J,length;
	double curForce = 0.0f;
	PLANE_4N_2D_Form *el;
	for (uint16 i=0; i < elems.size(); i++)
	{
		matX.zero();
		matS.zero();
		el = (PLANE_4N_2D_Form *) storage->el_form[elems[i]-1];
		matX[0][0] = 1+el->O[4][0];
		matX[0][1] = el->O[4][1];
		matX[1][0] = el->O[4][2];
		matX[1][1] = 1+el->O[4][3];

		matS[0][0] = el->S[4][0];
		matS[0][1] = el->S[4][2];
		matS[1][0] = el->S[4][2];
		matS[1][1] = el->S[4][1];
		
		J= matX[0][0]*matX[1][1] - matX[0][1]*matX[1][0];

		matT = matX * matS * matX.transpose() * (1.0f/J);
		double xx=storage->nodes[nodes1[i]-1].coord[0]+(*storage->vecQsum)[nodes1[i]*Node::nDOF-1] - 
					storage->nodes[nodes2[i]-1].coord[0]-(*storage->vecQsum)[nodes2[i]*Node::nDOF-1];
		double yy=storage->nodes[nodes1[i]-1].coord[1]+(*storage->vecQsum)[nodes1[i]*Node::nDOF-1+1] - 
					storage->nodes[nodes2[i]-1].coord[1]-(*storage->vecQsum)[nodes2[i]*Node::nDOF-1+1];
		length=sqrt(xx*xx+yy*yy);
		curForce += -matT[1][1]*length;
	}
	forces.push_back(curForce);

}

void Force_by_stress_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	echolog("***POST%d:%s OUTPUT:***\nLoadstep\tForce",nPost_proc,name.c_str());
	for (uint16 i=0; i < forces.size(); i++)
	{
		echolog("%d.\t%f",i,forces[i]);
	}
	echolog("***END POST%d:%s OUTPUT***",nPost_proc,name.c_str());

	//запись в файл
	char filename[100];
	sprintf_s(filename,100,"POST%d%s.txt",nPost_proc,name.c_str());
	ofstream file(filename,ios::trunc);
	if (!file)
	{
		warning("Stress_proc::post: Can't create a file with name %s", filename);
		return;
	}
	for (uint16 i=0; i < forces.size(); i++)
		file << forces[i] << endl; 
	file.close();
}



///////////////////////////////
Force_proc::Force_proc(FE_Storage *st) : Post_proc(st)
{
	name ="Force";
	active = true;
}

void Force_proc::pre(uint16 qLoadstep)
{
	/*echolog("STRESS_PROC INFORMATION:\nEl\tn1\tn2");
	for (uint16 i=0; i < elems.size(); i++)
	{
		echolog("%d.\t%d\t%d\t%d",i+1,elems[i],nodes1[i],nodes2[i]);
	}*/
	if (!active) return;
	forces.clear();
	forces.push_back(0.0f);
}

void Force_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	// пока реализуем только для ПНС
	if (storage->fem_type != FEM_PLANE && storage->fem_type != FEM_AXYM && storage->fem_type != FEM_PLANE_MIXED && storage->fem_type != FEM_PLANE_MIXED)
	{
		warning("Force_proc::process: this post_proc works only with PLANE and AXYM analysis");
		return;
	}

	storage->matK->zero();
	storage->vecF->zero();
	for (uint32 el = 0; el < storage->en; el++)
	{
		storage->el_form[el]->temp_k=0.0;
		storage->el_form[el]->build(el, storage);
		storage->el_form[el]->temp_k=1.0;
	}

	double curForce = 0.0f;
	for (uint16 i=0; i < nodes.size(); i++)
	{
		curForce+=(*storage->vecF)[nodes[i]*Node::nDOF-1+dofs[i]];
	}

	forces.push_back(curForce);
}

void Force_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	echolog("***POST%d:%s OUTPUT:***\nLoadstep\tForce",nPost_proc,name.c_str());
	for (uint16 i=0; i < forces.size(); i++)
	{
		echolog("%d.\t%f",i,forces[i]);
	}
	echolog("***END POST%d:%s OUTPUT***",nPost_proc, name.c_str());

	//запись в файл
	char filename[100];
	sprintf_s(filename,100,"POST%d %s.txt",nPost_proc,name.c_str());
	ofstream file(filename,ios::trunc);
	if (!file)
	{
		warning("Force_proc::post: Can't create a file with name %s", filename);
		return;
	}
	for (uint16 i=0; i < forces.size(); i++)
		file << forces[i] << endl; 
	file.close();
}

void Force_proc::load_from_file (const char *filename)
{
	ifstream file(filename);
	if (!file)
	{
		warning("Force_proc::load_from_file: Can't open input file `%s`",filename);
		return;
	}
	char buf[1024]="";
	uint32 val;
	uint16 key;
	while (file.getline(buf, 1024))
	{
		sscanf(buf, "%d %d", &key, &val);
		if (val < 1) break;
		if ((key < 0) && (key > 1)) break;
		nodes.push_back(val);
		dofs.push_back(key);
	}
	file.close();
}


////////////////////////////////////
RForce_proc::RForce_proc(FE_Storage *st) : Post_proc(st)
{
	name ="RForce";
	active = true;
}

void RForce_proc::pre(uint16 qLoadstep)
{
	if (!active) return;
	fx.clear();
	fy.clear();
}

void RForce_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	// пока реализуем только для ПНС
	if (curLoadstep != qLoadstep) return;

	storage->matK->zero();
	storage->vecF->zero();
	for (uint32 el = 0; el < storage->en; el++)
	{
		storage->el_form[el]->temp_k=0.0;
		storage->el_form[el]->build(el, storage);
		storage->el_form[el]->temp_k=1.0;
	}

	double curForce = 0.0f;
	for (uint16 i=0; i < nodes.size(); i++)
	{
		fx.push_back((*storage->vecF)[nodes[i]*Node::nDOF-1+0]);
		fy.push_back((*storage->vecF)[nodes[i]*Node::nDOF-1+1]);
	}
}


void RForce_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;

	//запись в файл
	char filename[100];
	sprintf_s(filename,100,"POST%d %s.txt",nPost_proc,name.c_str());
	ofstream file(filename,ios::trunc);
	if (!file)
	{
		warning("RForce_proc::post: Can't create a file with name %s", filename);
		return;
	}
	file << nodes.size() << endl;
	for (uint16 i=0; i < nodes.size(); i++)
		file << nodes[i] << "  ";
	file << endl;
	for (uint16 i=0; i < nodes.size(); i++)
		file << (storage->nodes)[nodes[i]-1].coord[0] << "  ";
	file << endl;
	for (uint16 i=0; i < nodes.size(); i++)
		file << (storage->nodes)[nodes[i]-1].coord[1] << "  ";
	file << endl;
	for (uint16 i=0; i < nodes.size(); i++)
		file << fx[i] << "  ";
	file << endl;
	for (uint16 i=0; i < nodes.size(); i++)
		file << fy[i] << "  ";
	file << endl;
	file.close();
}

void RForce_proc::load_from_file (const char *filename)
{
	ifstream file(filename);
	if (!file)
	{
		warning("RForce_proc::load_from_file: Can't open input file `%s`",filename);
		return;
	}
	char buf[1024]="";
	uint32 val;
	while (file.getline(buf, 1024))
	{
		sscanf(buf, "%d", &val);
		if (val < 1) break;
		nodes.push_back(val);
	}
	file.close();
}
#include "Strain_proc.h"

Strain_proc::Strain_proc(FE_Storage *st, string code) : Post_proc(st) 
{
	str_code = E_UNDEF;
	name ="Strain_proc";
	render = NULL;
	scene = NULL;
	char buf[200];
	if (code == "EX")
		str_code = E_X;
	else if (code == "EY")
		str_code = E_Y;
	else if (code == "EZ")
		str_code = E_Z;
	else if (code == "EXY")
		str_code = E_XY;
	else if (code == "EXZ")
		str_code = E_XZ;
	else if (code == "EYZ")
		str_code = E_YZ;
	else if (code == "ER")
		str_code = E_R;
	else if (code == "ET")
		str_code = E_T;
	else if (code == "ERZ")
		str_code = E_RZ;
	else if (code == "ERT")
		str_code = E_RT;
	else if (code == "EZT")
		str_code = E_ZT;
	else if (code == "EVOL")
		str_code = E_VOL;
	else {
		warning ("Unkown strain label %s", code.c_str());
		failed = true;
		return;
	}
	scene->label = "***" + code + "***";
	sprintf_s(buf, 200, "POST%d:%s(%s)",nPost_proc, name.c_str(),code.c_str());
	render = new Render(buf, 600,600);
	scene = new ColorE_Scene(); //TODO
	el_vals.clear();
	el_vals.assign(storage->en, 0.0f);
	render->attach(scene);
	//scene->make(st); //TODO
	scene->readMesh(st);
	read_Color(st);
	render->start();
	active = true;
}
Strain_proc::~Strain_proc()
{
	active = false;
	failed = true;
	if (render) render->stop();
	while (render->isStarted()) { };
	delete render;
	delete scene;
}

void Strain_proc::pre (uint16 qLoadstep)
{
	if (!active) return;
	el_vals.clear();
	el_vals.assign(storage->en, 0.0f);	
}
void Strain_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	WaitForSingleObject(scene->scene_lock,INFINITE);
	scene->readMesh(storage);
	read_Color(storage);
	ReleaseMutex(scene->scene_lock);
	
}
void Strain_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	//vector<uint16> rows;
	//vector<double> res;
	//rows.clear();
	//rows.push_back(1);
	//rows.push_back(15);
	//rows.push_back(50);
	//rows.push_back(75);
	//rows.push_back(100);
	//Mat<3,3> matX;
	//double val;
	//PLANE_4N_2D_Form *el_pl;
	//AXYM_4N_2D_Form *el_ax;
	//if (strain_code == E_VOL)
	//{
	//	char filename[100];
	//	sprintf_s(filename,100,"POST%d%s.txt",nPost_proc,name.c_str());
	//	ofstream file(filename,ios::trunc);
	//	if (!file)
	//	{
	//		warning("Stress_proc::post: Can't create a file with name %s", filename);
	//		return;
	//	}
	//	file << "Вывод объемных деформаци по сечениям после последнего шага решения" << endl;
	//	for (uint16 i=0; i < rows.size(); i++)
	//	{
	//		res.clear();
	//		uint16 row = rows[i];
	//		for (uint16 col=1; col < 41/*костыль*/; col++)
	//		{
	//			matX.zero();
	//			if (storage->fem_type == FEM_PLANE)
	//			{
	//				el_pl = (PLANE_4N_2D_Form *) storage->el_form[(row-1)*40+col-1];
	//				matX[0][0] = 1+2*el_pl->E[4][0];
	//				matX[0][1] = el_pl->E[4][2];
	//				matX[1][0] = el_pl->E[4][2];
	//				matX[1][1] = 1+2*el_pl->E[4][1];
	//				matX[2][2] = 1;
	//			}
	//			else if (storage->fem_type == FEM_AXYM)
	//			{
	//				el_ax = (AXYM_4N_2D_Form *) storage->el_form[(row-1)*40+col-1];
	//				matX[0][0] = 1+2*el_ax->E[4][0];
	//				matX[0][1] = el_ax->E[4][2];
	//				matX[1][0] = el_ax->E[4][2];
	//				matX[1][1] = 1+2*el_ax->E[4][1];
	//				matX[2][2] = 1+2*el_ax->E[4][3];
	//			}
	//			val= sqrt((matX[0][0]*matX[1][1] - matX[0][1]*matX[1][0])*matX[2][2]);
	//			res.push_back(val-1);
	//		}//цикл по cols
	//		file << "Сечение номер " << row << endl << "[";
	//		for (uint16 i=0; i < res.size(); i++)
	//			if (i!=res.size()-1)
	//				file << res[i] << ';';
	//			else
	//				file << res[i] << ']'<< endl;

	//		file << endl << endl << endl;
	//	}//цикл по сечениям
	//	file.close();
	//}//if statement
}

void Strain_proc::read_Color(FE_Storage *st)
{
	double val,min,max;
	if (storage->el_form.size() == 0 || str_code == E_UNDEF)
	{
		//дефолтный цвет
		scene->colors.assign(storage->en, -1);
		return;
	}
	//проверка: подходящая ли конфигурация
	if (_4N_2D_Form::nInt != 3)
	{
		warning("Strain_proc::process: this post_proc works only with integration scheme 3x3");
		return;
	}
	
	for (uint16 i=0; i < storage->en; i++)
	{
		el_vals[i] = storage->el_form[i]->getStrain(4,str_code);
		if (i==0) 
			min = val, max = val;
		else
		{
			if (val < min) min = val;
			if (val > max) max = val;
		}
	}
	scene->cmap.min = min;
	scene->cmap.max = max;
	for (uint16 i=0; i < storage->en; i++)
	{
		scene->colors[i] = scene->cmap.getIndex(el_vals[i]);
	}
}
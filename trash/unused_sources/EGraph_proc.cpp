#include "EGraph_proc.h"

EGraph_proc::EGraph_proc(FE_Storage *st, string str_code, uint16 gp) : Post_proc(st) 
{
	code = S_UNDEF;
	calc_scheme = gp;
	name ="EGraph_proc";
	char buf[200];
	scene = NULL;
	render = NULL;

	if (str_code == "EX")
		code = E_X;
	else if (str_code == "EY")
		code = E_Y;
	else if (str_code == "EZ")
		code = E_Z;
	else if (str_code == "EXY")
		code = E_XY;
	else if (str_code == "EXZ")
		code = E_XZ;
	else if (str_code == "EYZ")
		code = E_YZ;
	else if (str_code == "ER")
		code = E_R;
	else if (str_code == "ET")
		code = E_T;
	else if (str_code == "ERZ")
		code = E_RZ;
	else if (str_code == "ERT")
		code = E_RT;
	else if (str_code == "EZT")
		code = E_ZT;
	else if (str_code == "EVOL")
		code = E_VOL;
	else if (str_code == "SX")
		code = S_X;
	else if (str_code == "SY")
		code = S_Y;
	else if (str_code == "SZ")
		code = S_Z;
	else if (str_code == "SXY")
		code = S_XY;
	else if (str_code == "SXZ")
		code = S_XZ;
	else if (str_code == "SYZ")
		code = S_YZ;
	else if (str_code == "SR")
		code = S_R;
	else if (str_code == "ST")
		code = S_T;
	else if (str_code == "SRZ")
		code = S_RZ;
	else if (str_code == "SRT")
		code = S_RT;
	else if (str_code == "SZT")
		code = S_ZT;
	else if (str_code == "SP")
		code = S_P;
	else {
		warning ("Unkown label %s", str_code.c_str());
		failed = true;
		return;
	}
	//определим тип выводимой величины
	type = E_COMPONENT;
	name += "("+str_code+")"; 
	
	sprintf_s(buf, 200, "POST%d:%s",nPost_proc, name.c_str());
	scene = new ColorE_Scene(); //TODO
	scene->label = "***" + str_code + "***";
	render = new Render(buf, 600,600);
	el_vals.clear();
	el_vals.assign(storage->en, 0.0f);
	render->attach(scene);
	scene->readMesh(st);
	fillColors();
	render->start();
	active = true;
}
EGraph_proc::~EGraph_proc()
{
	active = false;
	failed = true;
	if (render) render->stop();
	echolog("start waiting of isStarted()");
	while (render->isStarted()) { };
	echolog("end waiting of isStarted()");
	delete render;
	delete scene;
}

void EGraph_proc::pre (uint16 qLoadstep)
{
	//if (!active) return;
	el_vals.clear();
	el_vals.assign(storage->en, 0.0f);
}
void EGraph_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	WaitForSingleObject(scene->scene_lock,INFINITE);
	scene->readMesh(storage);
	fillColors();
	ReleaseMutex(scene->scene_lock);
}
void EGraph_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
}

void EGraph_proc::fillColors ()
{
	double val,min,max;
	//в случае если еще задача не решалась
	if (storage->el_form.size() == 0)
	{
		//дефолтный цвет
		scene->colors.assign(storage->en, -1);
		return;
	}
	//проверка: подходящая ли конфигурация
	if (calc_scheme != GP_MEAN && calc_scheme >= _4N_2D_Form::nInt*_4N_2D_Form::nInt)
	{
		warning("EGraph_proc::fillColors: calc_scheme isn't appropriate for the current integration scheme");
		return;
	}
	
	for (uint16 i=0; i < storage->en; i++)
	{
		val=0.0f;
		uint16 g_start;
		uint16 g_end;
		if (calc_scheme == GP_MEAN)
		{
			g_start = 0;
			g_end = _4N_2D_Form::nInt*_4N_2D_Form::nInt;
		} else
		{
			g_start = calc_scheme;
			g_end = g_start+1;
		}
		for (uint16 j = g_start; j < g_end ; j++)
		{
			if (type == E_COMPONENT)
				val += storage->el_form[i]->getComponent(j,code, i, storage);
			else
			{
				warning("EGraph:fillColors: invalid type of element's output (%d)", type);
				return;
			}
		}
		el_vals[i] = val/(float)(g_end-g_start);
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

string EGraph_proc::getStatus ()
{
	return Post_proc::getStatus();
}
#include "Stress_proc.h"
//#include "ColorE_Scene.h"
Stress_proc::Stress_proc(FE_Storage *st, string code) : Post_proc(st) 
{
	str_code = S_UNDEF;
	name ="Stress_proc";
	char buf[200];
	scene = NULL;
	render = NULL;
	
	if (code == "SX")
		str_code = S_X;
	else if (code == "SY")
		str_code = S_Y;
	else if (code == "SZ")
		str_code = S_Z;
	else if (code == "SXY")
		str_code = S_XY;
	else if (code == "SXZ")
		str_code = S_XZ;
	else if (code == "SYZ")
		str_code = S_YZ;
	else if (code == "SR")
		str_code = S_R;
	else if (code == "ST")
		str_code = S_T;
	else if (code == "SRZ")
		str_code = S_RZ;
	else if (code == "SRT")
		str_code = S_RT;
	else if (code == "SZT")
		str_code = S_ZT;
	else {
		warning ("Unkown stress label %s", code.c_str());
		failed = true;
		return;
	}
	scene->label = "***" + code + "***";
	sprintf_s(buf, 200, "POST%d:%s(%s)",nPost_proc, name.c_str(),code.c_str());
	scene = new ColorE_Scene(); //TODO
	render = new Render(buf, 600,600);
	el_vals.clear();
	el_vals.assign(storage->en, 0.0f);
	render->attach(scene);
	scene->readMesh(st);
	read_Color(st);
	render->start();
	active = true;
}
Stress_proc::~Stress_proc()
{
	active = false;
	failed = true;
	if (render) render->stop();
	while (render->isStarted()) { };
	delete render;
	delete scene;
}

void Stress_proc::pre (uint16 qLoadstep)
{
	if (!active) return;
	el_vals.clear();
	el_vals.assign(storage->en, 0.0f);
}
void Stress_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	WaitForSingleObject(scene->scene_lock,INFINITE);
	scene->readMesh(storage);
	read_Color(storage);
	ReleaseMutex(scene->scene_lock);
}
void Stress_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
}

void Stress_proc::read_Color (FE_Storage *st)
{
	// код выводимого напряжения
	uint16 str_code=0;
	switch (stress_code)
	{
	case S_X:
		str_code = 1;
		break;
	case S_Y:
		str_code = 2;
		break;
	case S_XY:
		str_code = 3;
		break;
	case S_T:
		str_code = 4;
		break;
	default:
		warning("Stress_proc::process: unknown stress code!");
		return;
	}
	double val,min,max;
	//в случае если еще задача не решалась
	if (storage->el_form.size() == 0)
	{
		//дефолтный цвет
		scene->colors.assign(storage->en, -1);
		return;
	}
	//проверка: подходящая ли конфигурация
	if (_4N_2D_Form::nInt != 3)
	{
		warning("Stress_proc::process: this post_proc works only with integration scheme 3x3");
		return;
	}
	
	for (uint16 i=0; i < storage->en; i++)
	{
		el_vals[i] = storage->el_form[i]->getStress(4,str_code);
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
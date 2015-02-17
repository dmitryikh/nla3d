#include "Def_Form_proc.h"
#include "Def_Scene.h"
Def_Form_proc::Def_Form_proc(FE_Storage *st) : Post_proc(st) 
{
	name ="Def_form";
	char buf[200];
	sprintf_s(buf, 200, "POST%d:%s",nPost_proc, name.c_str());
	render = new Render(buf, 600,600);
	scene = new Def_Scene();
	render->attach(scene);
	scene->make(st);
	render->start();
	active = true;
}
Def_Form_proc::~Def_Form_proc()
{
	render->stop();
	echolog("start waiting of isStarted()");
	while (render->isStarted()) { Sleep(10);};
	echolog("end waiting of isStarted()");
	delete render;
	delete scene;
}

void Def_Form_proc::pre (uint16 qLoadstep)
{
}
void Def_Form_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (active)
		scene->make(storage);
}
void Def_Form_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{

}
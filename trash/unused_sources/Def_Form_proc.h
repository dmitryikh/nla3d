#pragma once
#include "sys.h"
#include "post_proc.h"
#include "Render.h"
#include "Scene.h"
#include "FE_Storage.h"

class Def_Form_proc : public Post_proc
{
public:
	Def_Form_proc(FE_Storage *st);
	virtual ~Def_Form_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
private:
	Render *render;
	Def_Scene *scene;
};
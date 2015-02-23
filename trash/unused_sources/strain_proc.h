#pragma once
#include "sys.h"
#include <vector>
#include "post_proc.h"
#include "Render.h"
#include "Scene.h"
#include "ColorE_Scene.h"
#include "FE_Storage.h"

//1-EX 2-EY 3-EZ 4-EXY 5-EXZ 6-EYZ
//7-ER 3-EZ 8-ET 9-ERZ 10-ERT 11-EZT
//12-EVOL

class Strain_proc : public Post_proc
{
public:
	Strain_proc(FE_Storage *st, string code);
	virtual ~Strain_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
private:
	uint16 str_code;
	string str_name;
	void read_Color (FE_Storage *st);
	Render *render;
	ColorE_Scene *scene;
	vector<float> el_vals;
};
#pragma once
#include "sys.h"
#include <vector>
#include "post_proc.h"
#include "Render.h"
#include "Scene.h"
#include "ColorE_Scene.h"
#include "FE_Storage.h"

//что выводит процессор
#define E_COMPONENT 1
//процессор графически выводит данные об элементе
class EGraph_proc : public Post_proc
{
public:
	EGraph_proc(FE_Storage *st, string code, uint16 gp);
	virtual ~EGraph_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
	string getStatus();
private:
	void fillColors ();
	uint16 code;
	uint16 type;
	uint16 calc_scheme;
	Render *render;
	ColorE_Scene *scene;
	vector<float> el_vals;
};
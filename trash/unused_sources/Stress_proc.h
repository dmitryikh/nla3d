#pragma once
#include "sys.h"
#include <vector>
#include "post_proc.h"
#include "Render.h"
#include "Scene.h"
#include "ColorE_Scene.h"
#include "FE_Storage.h"

#define GP_MEAN 100 //среднее значение по элементу
//что выводит процессор
#define E_STRESS 1
#define E_STRAIN 2
//процессор графически выводит данные об элементе
class EGraph_proc : public Post_proc
{
public:
	EGraph_proc(FE_Storage *st, string code, uint16 gp);
	virtual ~EGraph_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
	virtual string getStatus();
private:
	void read_Color (FE_Storage *st);
	uint16 code;
	uint16 type;
	Render *render;
	ColorE_Scene *scene;
	vector<float> el_vals;
};
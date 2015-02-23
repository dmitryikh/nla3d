#pragma once
#include "sys.h"
#include "Colormap.h"

extern Color el_line;
extern Color el_fill;
extern Color bc_zero;
extern Color bc_non_zero;
extern Color background;

class Scene 
{
public:
	virtual void draw ()=0;
	virtual void setup_view (uint16 width, uint16 height)=0;
	uint16 font_base;
private:
};

struct Scene_ubound
{
public:
	float fi;
	double x;
	double y;
	float len;
	Color color;
};


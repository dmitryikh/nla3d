#pragma once
#include "sys.h"
#include <gl\GL.h>
#include <vector>
#include "Scene.h"
#include "Colormap.h"
#include "math\Vec.h"
#include "FE_Storage.h"


class ColorN_Scene : public Scene
{
public:
	ColorN_Scene() {
		scene_lock=CreateMutex(NULL, FALSE, NULL);
		wide=1.0f;
	}
	HANDLE scene_lock;
	vector<Vec<2>> points[4];
	vector<int> colors;
	double wide;
	Colormap cmap;
	void readMesh (FE_Storage *st);
	virtual void draw ();
	virtual void setup_view (uint16 width, uint16 height);
	uint16 width;
	uint16 height;
	string label;
};
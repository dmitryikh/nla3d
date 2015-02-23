#pragma once
#include "sys.h"
#include <gl\GL.h>
#include "Scene.h"
#include "FE_Storage.h"
#include "math\Vec.h"
class FE_Scene : public Scene {
public:
	FE_Scene () {
		scene_lock=CreateMutex(NULL, FALSE, NULL);
		wide=1.0f;
	}
	HANDLE scene_lock;
	vector<Vec<2>> p1;
	vector<Vec<2>> p2;
	vector<Vec<2>> quads;
	vector<Scene_ubound> bounds;
	double wide;
	void make (FE_Storage *storage);
	virtual void draw ();
	virtual void setup_view (uint16 width, uint16 height);
};
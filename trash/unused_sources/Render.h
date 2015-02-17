#pragma once
#include "sys.h"
#include <vector>
#include <gl/GL.h>
#include <gl/GLU.h>
#include "math\Vec.h"
#include "Scene.h"

class Render
{
public:
	Render (const char* wnd_name, uint16 w, uint16 h);
	void attach(Scene* s) {
		assert(s);
		scene=s;
	}
	static LRESULT WINAPI WndProc (HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
	void loop ();
	inline void start () {
		_beginthread(fork,0,this);
	}
	inline void stop () {
		close_signal=true;
	}
	bool isStarted() {
		return started;
	};

	Vec<2> trans;
	double scale;
	~Render();
private:
	static void fork (void *ptr);
	void changeSize (uint16 w, uint16 h);
	void render_scene();
	void setDCPixelFormat();
	const uint16 this_render;
	static HANDLE render_lock; //назначение мьютекса - чтобы в одно время поток openGL занимал только один поток программы
	static uint16 render_count;
	Scene *scene;
	uint32 width, height;
	string name;
	bool started;
	bool close_signal;
	HDC hDC;
	HGLRC hRC;

	GLuint font_base;
};

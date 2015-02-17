#include "Render.h"

HANDLE Render::render_lock = CreateMutex(NULL, FALSE, NULL);
uint16 Render::render_count=0; // количество запущенных рендеров 
Render::Render(const char* wnd_name, uint16 w, uint16 h):this_render(render_count++) {
	width=w;
	height=h;
	scene=NULL;
	name.assign(wnd_name);
	hDC=NULL;
	hRC=NULL;
	started=false;
	close_signal=false;
	scale=1.0f;
	font_base = 0;
}

Render::~Render () {
	echolog("~Render");
	render_count--;
}

LRESULT WINAPI Render::WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	if(message == WM_NCCREATE){
		SetWindowLong(hWnd, GWL_USERDATA, (LONG)((LPCREATESTRUCT)lParam)->lpCreateParams);
	}

	Render *render=NULL;
	render = (Render*)GetWindowLong(hWnd, GWL_USERDATA);
	/*if (!render)
		warning("Can't obtain render pointer in static WndProc");*/
	//обработка!
	//return pWindow->WndProc(hwnd, uMsg, wParam, lParam);
	switch (message)
	{
	case  WM_CREATE:
		render->hDC = GetDC(hWnd);
		render->setDCPixelFormat();
		WaitForSingleObject(render_lock, INFINITE );
		render->hRC=wglCreateContext(render->hDC);
		wglMakeCurrent(render->hDC, render->hRC);
		ReleaseMutex(render_lock);
		SetTimer(hWnd, 1, 33, NULL);
		break;
	case	WM_DESTROY:
		WaitForSingleObject(render_lock, INFINITE );
		wglDeleteContext(render->hRC);
		ReleaseMutex(render_lock);
		PostQuitMessage(0);
		break;
	case  WM_SIZE:
		render->changeSize(LOWORD(lParam), HIWORD(lParam));
		break;
	case WM_TIMER:
		InvalidateRect(hWnd,NULL,FALSE);
		break;
	case WM_PAINT:	
		render->render_scene();
		SwapBuffers(render->hDC);
		ValidateRect(hWnd, NULL);
		break;
	case WM_KEYDOWN:
		switch (wParam) {
		case VK_LEFT:
			render->trans[0]-=-1.0f*render->scale;
			break;
		case VK_RIGHT:
			render->trans[0]+=-1.0f*render->scale;
			break;
		case VK_UP:
			render->trans[1]+=-1.0f*render->scale;
			break;
		case VK_DOWN:
			render->trans[1]-=-1.0f*render->scale;
			break;
		case VK_F1:
			render->scale*=1.1f;
			break;
		case VK_F2:
			render->scale/=1.1f;
			break;
		default:;
		}
		break;
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return (0L);
}

void Render::loop() 
{
	//создание окна и заупск цикла
	HINSTANCE proc = GetModuleHandle(NULL) ;
	WNDCLASS wc;
	HWND hwnd;
	MSG msg;
	wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	wc.lpfnWndProc = (WNDPROC) Render::WndProc;
	wc.cbClsExtra=0;
	wc.cbWndExtra=0;
	wc.hInstance=proc;
	wc.hIcon=NULL;
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground=NULL;
	wc.lpszMenuName=NULL;
	string str = "Render";	
	str.push_back('1'+this_render);
	wstring str1(str.begin(), str.end());
	wc.lpszClassName = str1.c_str();
	if (RegisterClass(&wc)==0)
		error("Can`t register window class %s",name.c_str());
	hwnd=CreateWindow((LPCWSTR) str1.c_str(),
		wstring(name.begin(), name.end()).c_str(),
		WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN | WS_CLIPSIBLINGS,
		0,0,
		width, height,
		NULL, NULL,
		proc,
		(LPVOID) this);
	if (hwnd==NULL)
		error("Can`t create window %s", name.c_str());
	debug("Launch %d render `%s`",this_render, name.c_str());
	ShowWindow(hwnd, SW_SHOW);
	UpdateWindow(hwnd);
	started=true;
	while (GetMessage(&msg, NULL,0,0))
	{
		//echolog("%s MSG!",name.c_str());
		if (close_signal==true)
		{
			echolog("kill signal No %s",name.c_str());
			KillTimer(hwnd, 1);
			SendMessage(hwnd, WM_DESTROY,0,0);
		}
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	echolog("started to false");
	started = false;
}

void Render::setDCPixelFormat()
{
	int nPixelFormat;

    static PIXELFORMATDESCRIPTOR pfd = {
        sizeof(PIXELFORMATDESCRIPTOR), // Size of this structure
        1,                             // Version of this structure
        PFD_DRAW_TO_WINDOW |           // Draw to window (not to bitmap)
        PFD_SUPPORT_OPENGL |           // Support OpenGL calls in window
        PFD_DOUBLEBUFFER,              // Double-buffered mode
        PFD_TYPE_RGBA,                 // RGBA color mode
        32,                            // Want 32-bit color
        0,0,0,0,0,0,                   // Not used to select mode
        0,0,                           // Not used to select mode
        0,0,0,0,0,                     // Not used to select mode
        16,                            // Size of depth buffer
        0,                             // Not used here
        0,                             // Not used here
        0,                             // Not used here
        0,                             // Not used here
        0,0,0 };                       // Not used here


    // Choose a pixel format that best matches that described in pfd
    nPixelFormat = ChoosePixelFormat(hDC, &pfd);

    // Set the pixel format for the device context
    SetPixelFormat(hDC, nPixelFormat, &pfd);
}

void Render::changeSize(uint16 w, uint16 h)
{
	//if (!started) return;
	width=w;
	height=h;
	//echolog("%s -> WAIT renderlock (change size)", name.c_str());
	WaitForSingleObject( render_lock, INFINITE );
	//echolog("%s -> SET renderlock (change size)", name.c_str());
	wglMakeCurrent(hDC, hRC);

	if (!font_base)
	{
		HFONT	font;										// Windows Font ID
		HFONT	oldfont;									// Used For Good House Keeping

		font_base = glGenLists(96);								// Storage For 96 Characters
		font = CreateFont(	-16,							// Height Of Font
							0,								// Width Of Font
							0,								// Angle Of Escapement
							0,								// Orientation Angle
							FW_BOLD,						// Font Weight
							FALSE,							// Italic
							FALSE,							// Underline
							FALSE,							// Strikeout
							ANSI_CHARSET,					// Character Set Identifier
							OUT_TT_PRECIS,					// Output Precision
							CLIP_DEFAULT_PRECIS,			// Clipping Precision
							ANTIALIASED_QUALITY,			// Output Quality
							FF_DONTCARE|DEFAULT_PITCH,		// Family And Pitch
							L"Arial");					// Font Name

		oldfont = (HFONT)SelectObject(hDC, font);           // Selects The Font We Want
		wglUseFontBitmaps(hDC, 32, 96, font_base);				// Builds 96 Characters Starting At Character 32
		SelectObject(hDC, oldfont);							// Selects The Font We Want
		DeleteObject(font);		
		if (scene) 
			scene->font_base = font_base;
	}

	glViewport(0,0,w,h);
	if (scene) scene->setup_view(w,h);
	wglMakeCurrent(NULL, NULL);
	ReleaseMutex(render_lock);
	//echolog("%s -> RELEASE renderlock (change size)", name.c_str());
}

void Render::render_scene()
{
	if (!started) return;
	//echolog("%s -> WAIT renderlock", name.c_str());
	WaitForSingleObject(render_lock, INFINITE );
	//echolog("%s -> SET renderlock", name.c_str());
	wglMakeCurrent(hDC, hRC);
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScaled(scale, scale,0);
	glTranslated(trans[0], trans[1],0.0f);
	//
	if (scene) scene->draw();
	wglMakeCurrent(NULL, NULL);
	ReleaseMutex(render_lock);
	//echolog("%s -> RELEASE renderlock", name.c_str());
}

void Render::fork (void *ptr)
{
	((Render*)ptr)->loop();
}
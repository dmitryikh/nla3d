#include "ColorE_Scene.h"

void ColorE_Scene::readMesh (FE_Storage *st)
{
	double multiple = 1.0f; //множитель перемещений
	WaitForSingleObject( scene_lock, INFINITE );
	points[0].clear();
	points[1].clear();
	points[2].clear();
	points[3].clear();
	for (uint32 el=0; el < st->en; el++)
	{
		if (el==0)
			wide=fabs(st->nodes[st->elements[el].nodes[0]-1].coord[0]);
		for (uint32 i=0; i < 4; i++)
		{
			Vec<2> coord;
			coord[0] = st->nodes[st->elements[el].nodes[i]-1].coord[0];
			coord[1] = st->nodes[st->elements[el].nodes[i]-1].coord[1];
			if (fabs(coord[0]) > wide) wide=fabs(coord[0]);
			if (fabs(coord[1]) > wide) wide=fabs(coord[1]);
			if (st->vecQsum)
			{
				coord[0]+=(*st->vecQsum)[st->elements[el].nodes[i]*Node::nDOF-1] * multiple;
				coord[1]+=(*st->vecQsum)[st->elements[el].nodes[i]*Node::nDOF-1+1] * multiple;
			}
			points[i].push_back(coord);
		}
	}
	//сначала всё - серый цвет
	colors.assign(st->en,-1);
	ReleaseMutex(scene_lock);
}

void ColorE_Scene::draw () 
{
	WaitForSingleObject( scene_lock, INFINITE ); //TIP: нужно лочить класс, чтобы никто не зименил данные во время прорисовки кадра
	glLineWidth(1.0f);
	glDisable(GL_DEPTH_TEST);
	glBegin(GL_QUADS);
	for (uint32 i=0; i<points[0].size(); i++)
	{
		glColor3fv(cmap[colors[i]].rgb);
		for (uint32 j=0; j < 4; j++)
			glVertex2d(points[j][i][0],points[j][i][1]);
	}
	glEnd();
	glColor3fv(el_line.rgb);
	glBegin(GL_LINES);
	for (uint32 i=0; i<points[0].size(); i++)
	{
		glVertex2d(points[0][i][0], points[0][i][1]);
		glVertex2d(points[1][i][0], points[1][i][1]);

		glVertex2d(points[1][i][0], points[1][i][1]);
		glVertex2d(points[2][i][0], points[2][i][1]);

		glVertex2d(points[2][i][0], points[2][i][1]);
		glVertex2d(points[3][i][0], points[3][i][1]);

		glVertex2d(points[3][i][0], points[3][i][1]);
		glVertex2d(points[0][i][0], points[0][i][1]);
	}
	glEnd();
	glLineWidth(2.0f);
	glLineStipple(5, 0xF6F6);
	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex2f(0.0f, 0.0f);
	glVertex2f(5.0f, 0.0f);
	glEnd();
	glEnable(GL_LINE_STIPPLE);
	
	glBegin(GL_LINES);
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex2f(0.0f, -10.0f);
	glVertex2f(0.0f, 10.0f);
	glEnd();
	glDisable(GL_LINE_STIPPLE);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	double ratio = (float)height/(float)width;
	glOrtho(0.0f, width,0.0f,height,1.0f,-1.0f);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	double leftborder  = 40.0f;
	double bottomborder = 30.0f;
	double length = width-2*leftborder;
	double up = 20.0f;
	glLineWidth(1.0f);
	char str[100];
	glPushAttrib(GL_LIST_BIT);							// Pushes The Display List Bits
		glListBase(this->font_base- 32);	// Sets The Base Character to 32
	for(uint16 i=0; i < this->cmap.colors.size(); i++)
	{
		glBegin(GL_QUADS);
		glColor3fv(cmap.colors[i].rgb);
		glVertex2d(leftborder+length/cmap.colors.size()*i, bottomborder);
		glVertex2d(leftborder+length/cmap.colors.size()*i, bottomborder+up);
		glVertex2d(leftborder+length/cmap.colors.size()*(i+1), bottomborder+up);
		glVertex2d(leftborder+length/cmap.colors.size()*(i+1), bottomborder);
		glEnd();
		glColor3d(0.0f,0.0f,0.0f);
		glBegin(GL_LINE_LOOP);
		glVertex2d(leftborder+length/cmap.colors.size()*i, bottomborder);
		glVertex2d(leftborder+length/cmap.colors.size()*i, bottomborder+up);
		glVertex2d(leftborder+length/cmap.colors.size()*(i+1), bottomborder+up);
		glVertex2d(leftborder+length/cmap.colors.size()*(i+1), bottomborder);
		glEnd();
		sprintf_s(str, 100, "%7.3f",cmap.min + i*(cmap.max-cmap.min)/cmap.colors.size());
		glRasterPos2f(leftborder+length/cmap.colors.size()*i-20, bottomborder-16);							
		glCallLists(strlen(str), GL_UNSIGNED_BYTE, str);	// Draws The Display List Text
			
	}

	sprintf_s(str, 100, "%7.3f",cmap.max);
	glRasterPos2f(leftborder+length-20, bottomborder-16);							
	glCallLists(strlen(str), GL_UNSIGNED_BYTE, str);	// Draws The Display List Text
	glRasterPos2f(width/2-20, height-50);
	glCallLists(label.length(), GL_UNSIGNED_BYTE, label.c_str());
	glPopAttrib();

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glFlush();

	ReleaseMutex( scene_lock );
}

void ColorE_Scene::setup_view (uint16 width, uint16 height)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	this->width = width;
	this->height = height;
	double plus=wide*0.1f;
	double ratio = (double)width/(double)height;
	if (width < height)
		glOrtho(-wide-plus, wide+plus,-(wide+plus)/ratio,(wide+plus)/ratio,1.0f,-1.0f);
	else
		glOrtho(-(wide+plus)*ratio, (wide+plus)*ratio,-(wide+plus),(wide+plus),1.0f,-1.0f);
	glMatrixMode(GL_MODELVIEW);
	glClearColor(background.rgb[0],background.rgb[1],background.rgb[2],0.0f);
	
}
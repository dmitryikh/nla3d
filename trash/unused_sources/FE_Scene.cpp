#include "FE_Scene.h"
#include "FE_Storage.h"

void FE_Scene::draw () 
{
	WaitForSingleObject( scene_lock, INFINITE );
	glColor3d(20/255.0f,220/255.0f,220/255.0f);
	glLineWidth(2.0f);
	glBegin(GL_QUADS);
	for (uint32 i=0; i<quads.size(); i++)
	{
		glVertex2d(quads[i][0], quads[i][1]);
	}
	glEnd();
	glColor3d(5/255.0f,180/255.0f,90/255.0f);
	glBegin(GL_LINES);
	for (uint32 i=0; i<p1.size(); i++)
	{
		glVertex2d(p1[i][0], p1[i][1]);
		glVertex2d(p2[i][0], p2[i][1]);
	}
	glEnd();
	double fi = 3.141593f/2.6f;
	float len=2.0f;
	const double pi = 3.141592;
	
	for (uint32 i=0; i<bounds.size(); i++)
	{
		glColor3fv(bounds[i].color.rgb);
		glBegin(GL_LINE_LOOP);
		glVertex2d(bounds[i].x, bounds[i].y);
		glVertex2d(bounds[i].x+bounds[i].len*cos(-(fi+bounds[i].fi)), bounds[i].y+bounds[i].len*sin(-(fi+bounds[i].fi)));
		glVertex2d(bounds[i].x+bounds[i].len*cos(-(pi-fi+bounds[i].fi)), bounds[i].y+bounds[i].len*sin(-(pi-fi+bounds[i].fi)));
		glEnd();
	}
	

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
	glFlush();
	ReleaseMutex( scene_lock );
}

void FE_Scene::setup_view (uint16 width, uint16 height)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	double plus=wide*0.1f;
	double ratio = (double)width/(double)height;
	if (width < height)
		glOrtho(-wide-plus, wide+plus,-(wide+plus)/ratio,(wide+plus)/ratio,1.0f,-1.0f);
	else
		glOrtho(-(wide+plus)*ratio, (wide+plus)*ratio,-(wide+plus),(wide+plus),1.0f,-1.0f);
	glMatrixMode(GL_MODELVIEW);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
}

void FE_Scene::make (FE_Storage *storage)
{
	WaitForSingleObject( scene_lock, INFINITE );
		p1.clear();
		p2.clear();
		quads.clear();
		for (uint32 el=0; el < storage->en; el++)
		{
			if (el==0)
				wide=fabs(storage->nodes[storage->elements[el].nodes[0]-1].coord[0]);
			for (uint32 i=0; i < 4; i++)
			{
				double x1,y1,x2,y2;
				x1 = storage->nodes[storage->elements[el].nodes[i]-1].coord[0];
				y1 = storage->nodes[storage->elements[el].nodes[i]-1].coord[1];

				if (i==3)
				{
					x2 = storage->nodes[storage->elements[el].nodes[0]-1].coord[0];
					y2 = storage->nodes[storage->elements[el].nodes[0]-1].coord[1];
				} 
				else 
				{
					x2 = storage->nodes[storage->elements[el].nodes[i+1]-1].coord[0];
					y2 = storage->nodes[storage->elements[el].nodes[i+1]-1].coord[1];
				}
				if (fabs(x1) > wide) wide=fabs(x1);
				if (fabs(x2) > wide) wide=fabs(x2);
				if (fabs(y1) > wide) wide=fabs(y1);
				if (fabs(y2) > wide) wide=fabs(y2);
				p1.push_back(Vec<2>(x1,y1));
				p2.push_back(Vec<2>(x2,y2));
			}
			quads.push_back(storage->nodes[storage->elements[el].nodes[0]-1].coord);
			quads.push_back(storage->nodes[storage->elements[el].nodes[1]-1].coord);
			quads.push_back(storage->nodes[storage->elements[el].nodes[2]-1].coord);
			quads.push_back(storage->nodes[storage->elements[el].nodes[3]-1].coord);	
		}
		for (uint32 i = 0; i < storage->bounds.size(); i++)
		{
			Scene_ubound sb;
			sb.fi=0;
			sb.x=storage->nodes[storage->bounds[i].node-1].coord[0];
			sb.y=storage->nodes[storage->bounds[i].node-1].coord[1];
			sb.len=(float) wide/25.0f;
			if (storage->bounds[i].key==D_UX) sb.fi= (float) 3.141593/2.0f;
			
			sb.color=bc_zero;
			if (fabs(storage->bounds[i].value) > 1.0e-7)
			{
				sb.color = bc_non_zero;
				if (storage->bounds[i].value < 0.0f) sb.fi+=(float) 3.141593f;
			}
			bounds.push_back(sb);
		}
	ReleaseMutex(scene_lock);
}

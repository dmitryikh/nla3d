#pragma once
#include "sys.h"
#include <vector>

struct Color
{
	Color() {
		rgb[0]=0.0f;
		rgb[1]=0.0f;
		rgb[2]=0.0f;
	}
	Color (float r, float g, float b) {
		rgb[0]=r;
		rgb[1]=g;
		rgb[2]=b;
	}
	Color& operator= (const Color &op)
	{
		this->rgb[0]=op.rgb[0];
		this->rgb[1]=op.rgb[1];
		this->rgb[2]=op.rgb[2];
		return *this;
	}
	float rgb[3];
};

const Color gray(0.5f, 0.5f, 0.5f);

class Colormap
{
public:
	Colormap ()
	{
		loadDefault();
		min = 0.0f;
		max = 0.0f;
	}
	vector<Color> colors;
	float min;
	float max;
	Color getColor (double val)
	{
		uint16 n = colors.size();
		if (abs(min-max) < 10e-8)
			return gray;
		int i = floor((val-min)*n/(max-min));
		if (i<0 || i > n-1)
			return gray;
		return colors[i];
	}
	int getIndex (double val)
	{
		uint16 n = colors.size();
		if (abs(min-max) < 10e-8)
			return -1;
		int i = floor((val-min)*n/(max-min));
		
		if (i==n) 
			i=n-1;
		if (i<0 || i > n-1)
			return -1;
		return i;
	}
	Color operator[] (int index)
	{
		assert(index < (int) colors.size());
		if (index < 0)
			return gray;
		return colors[index];
	}
	void loadDefault ()
	{
		colors.clear();
		colors.push_back(Color(0.0f, 0.0f, 1.0f)); //синий
		colors.push_back(Color(0.0f, 0.5f, 1.0f)); //голубой
		colors.push_back(Color(0.0f, 1.0f, 1.0f)); //салатовый
		colors.push_back(Color(0.5f, 1.0f, 0.5f)); //зеленый
		colors.push_back(Color(1.0f, 1.0f, 0.0f)); //желтый
		colors.push_back(Color(1.0f, 0.5f, 0.0f)); //оранжевый
		colors.push_back(Color(1.0f, 0.25f, 0.0f)); //оранж+красный
		colors.push_back(Color(1.0f, 0.0f, 0.0f)); //красный
	}
};
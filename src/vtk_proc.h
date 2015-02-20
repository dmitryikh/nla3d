#pragma once
#include "sys.h"
#include "post_proc.h"
#include "FE_Storage.h"

class Vtk_proc : public Post_proc
{
public:
	Vtk_proc(FE_Storage_Interface *st, string _fileName);
	virtual ~Vtk_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
private:
	string file_name;
	uint16 a_dim; // размерность задачи
	uint16 nInt; //схема интегрирования в задачи
	vector<el_component> comp_codes;
	void write_header (ofstream &file);
	void write_geometry(ofstream &file, bool def);
	void write_point_data(ofstream &file, bool zero);
	void write_pointed_el_component(ofstream &file, el_component code);
	void write_pointed_el_component2(ofstream &file, el_component code);
	void write_pointed_el_component3(ofstream &file, el_component code);
	void write_cell_data(ofstream &file, bool zero);
};

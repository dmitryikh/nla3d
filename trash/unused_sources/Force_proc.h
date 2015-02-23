#pragma once
#include "sys.h"
#include <vector>
#include "post_proc.h"
#include "FE_Storage.h"

class Force_by_stress_proc : public Post_proc
{
public:
	Force_by_stress_proc(FE_Storage *st);
	virtual ~Force_by_stress_proc() { };
	//~Stress_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
	vector<uint32> elems;
	vector<uint32> nodes1;
	vector<uint32> nodes2;
	vector<double> forces;
};

class Force_proc : public Post_proc
{
public:
	Force_proc(FE_Storage *st);
	virtual ~Force_proc() { };
	//~Stress_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
	void load_from_file (const char *filename);
	vector<uint32> nodes;
	vector<uint32> dofs;
	vector<double> forces;
};

class RForce_proc : public Post_proc
{
public:
	RForce_proc(FE_Storage *st);
	virtual ~RForce_proc() { };
	//~Stress_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
	void load_from_file (const char *filename);
	vector<uint32> nodes;
	vector<double> fx;
	vector<double> fy;
};



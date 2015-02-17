#pragma once
#include "sys.h"
#include "post_proc.h"
#include "FE_Storage.h"

class EList_proc : public Post_proc
{
public:
	EList_proc(FE_Storage *st, string str_code, uint16 gp);
	virtual ~EList_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
	string getStatus();
	vector<uint32> elems;
	uint16 every_n_ls;
private:
	uint16 calc_scheme;
	uint16 code;
	string file_name;
};
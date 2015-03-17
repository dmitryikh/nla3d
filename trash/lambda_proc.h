
#pragma once
#include "sys.h"
#include "post_proc.h"
#include "FE_Storage.h"
#include "element_MIXED_8N_3D_P0.h"

class lambda_proc : public Post_proc
{
public:
	lambda_proc(FE_Storage_Interface *st, string _fileName);
	virtual ~lambda_proc();
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
private:
	string file_name;
};

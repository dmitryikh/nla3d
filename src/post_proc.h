#pragma once
#include <string>
#include "sys.h"

class Def_Scene;
class FE_Storage_Interface;
//Data_Processor
class Post_proc 
{
public:
	Post_proc(FE_Storage_Interface *st);
	~Post_proc() { };
	virtual void pre (uint16 qLoadstep)=0;
	virtual void process (uint16 curLoadstep, uint16 qLoadstep)=0;
	virtual void post (uint16 curLoadstep, uint16 qLoadstep)=0;
	string getStatus ();
	uint16 getnPost_num () {
		return nPost_proc;
	}
	void setActive (bool act)
	{
		if (failed)
		{
			warning("post_proc %d (%s): can't set active, because post_proc has been already failed", nPost_proc, name.c_str());
			return;
		}
		active = act;
	}
	bool getActive ()
	{
		return active;
	}

	uint16 getProc_num ()
	{
		return nPost_proc;
	}
	template <typename el_type> friend class FE_Storage;
protected:
	FE_Storage_Interface *storage;
	uint16 nPost_proc;
	string name;
	bool active;
	bool failed;
};


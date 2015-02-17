#include "post_proc.h"
#include "FE_Storage.h"
Post_proc::Post_proc(FE_Storage_Interface *st) 
{
	assert(st);
	storage = st;
	storage->add_post_proc(this);
	active = false;
	failed = false;
}

string Post_proc::getStatus ()
{
	char buf[300];
	sprintf_s(buf, 200, "Post_proc No %d: %s\n\tActive - %s, Failed - %s",nPost_proc, name.c_str(),active?"true":"false",active?"true":"false");
	return string(buf);
}
#include "sys.h"
#include "post_proc.h"
#include "FE_Storage.h"


class Reaction_proc : public Post_proc
{
public:
	Reaction_proc(FE_Storage_Interface *st);
	virtual ~Reaction_proc() { };
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
	//void load_from_file (const char *filename);
	vector<uint32> nodes;
	vector<uint16> dofs;
	string filename;
};
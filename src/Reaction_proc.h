#include "sys.h"
#include "post_proc.h"
#include "FE_Storage.h"


class Reaction_proc : public Post_proc
{
public:
	Reaction_proc(FE_Storage *st);
	Reaction_proc(FE_Storage *st, string _filename);
	virtual ~Reaction_proc() { };
	virtual void pre (uint16 qLoadstep);
	virtual void process (uint16 curLoadstep, uint16 qLoadstep);
	virtual void post (uint16 curLoadstep, uint16 qLoadstep);
  vector<double> getReactions ();
	
	vector<uint32> nodes;
	vector<uint16> dofs;
protected:
  string filename;
	vector<double> reactVec;
};

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "FEStorage.h"

namespace nla3d {

class FEStorage;
class FESolver;

//Data_Processor
class PostProcessor {
public:
	PostProcessor(FEStorage *st);
	virtual ~PostProcessor() { };
	virtual void pre ()=0;
	virtual void process (uint16 curLoadstep)=0;
	virtual void post (uint16 curLoadstep)=0;
	std::string getStatus ();
	uint16 getnPost_num () {
		return nPost_proc;
	}
	void setActive (bool act) {
		if (failed) {
			LOG(WARNING) << "PostProcessor " << nPost_proc << " (" << name << "): can't set active,"
          << "because post_proc has been already failed";
			return;
		}
		active = act;
	}
	bool getActive () {
		return active;
	}

	uint16 getProc_num () {
		return nPost_proc;
	}

	friend class FEStorage;
	friend class FESolver;
protected:
	FEStorage *storage;
	uint16 nPost_proc;
	std::string name;
	bool active;
	bool failed;
};

} // namespace nla3d

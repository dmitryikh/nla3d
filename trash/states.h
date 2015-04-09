// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"

namespace nla3d {

class States {
public:
	States ();
	
	enum t_UI16 {
		UI16_CURINTPOINT = 0,
		UI16_EQUILITER,
		UI16_LASTID
	};
	
	enum t_UI32 {
		UI32_CURELEM = 0,
		UI32_CURLOADSTEP,
		UI32_CURSUBSTEP,
		UI32_LASTID
	};
	
	enum t_DOUBLE {
		DOUBLE_MATM = 0,
		DOUBLE_HYDPRES,
		DOUBLE_LASTID
	};
	enum t_PTR {
		PTR_CURMATER = 0,
		PTR_LASTID
	};

	public:
	double 	getdouble	(t_DOUBLE id);
	uint32 	getuint32	(t_UI32 id);
	uint16 	getuint16	(t_UI16 id);
	void* 	getptr		(t_PTR id);
	
	void 	setdouble	(t_DOUBLE id, double);
	void 	setuint32	(t_UI32 id, uint32);
	void 	setuint16	(t_UI16 id, uint16);
	void 	setptr		(t_PTR, void*);
	
	bool 	ifdefdouble	(t_DOUBLE id);
	bool 	ifdefuint32	(t_UI32 id);
	bool 	ifdefuint16	(t_UI16 id);
	bool 	ifdefptr 	(t_PTR id);
	
	void 	undefinedouble	(t_DOUBLE id);
	void 	undefineuint32	(t_UI32 id);
	void 	undefineuint16	(t_UI16 id);
	void 	undefineptr 	(t_PTR id);
	private:
	uint16 states_UI16 		[UI16_LASTID];
	uint32 states_UI32 		[UI32_LASTID];
	double states_DOUBLE 	[DOUBLE_LASTID];
	void *states_PTR		[PTR_LASTID];
	
	bool if_def_UI16		[UI16_LASTID];
	bool if_def_UI32		[UI32_LASTID];
	bool if_def_DOUBLE		[DOUBLE_LASTID];
	bool if_def_PTR			[PTR_LASTID];
	
};


extern States GlobStates;

} // namespace nla3d

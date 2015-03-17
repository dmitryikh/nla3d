// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "States.h"

namespace nla3d {

States GlobStates;

States::States () {
	uint16 i;
	for (i=0; i < UI16_LASTID; i++) {
			if_def_UI16[i] = false;
	}
	for (i=0; i < UI32_LASTID; i++) {
			if_def_UI32[i] = false;
	}
	for (i=0; i < DOUBLE_LASTID; i++) {
			if_def_DOUBLE[i] = false;
	}
	for (i=0; i < PTR_LASTID; i++) {
			if_def_PTR[i] = false;
	}
}

double States::getdouble(t_DOUBLE id) {
	assert(id < DOUBLE_LASTID);
	assert(if_def_DOUBLE[id]);
	return states_DOUBLE[id];
}

uint16 States::getuint16(t_UI16 id) {
	assert(id < UI16_LASTID);
	assert(if_def_UI16[id]);
	return states_UI16[id];
}

uint32 States::getuint32(t_UI32 id) {
	assert(id < UI32_LASTID);
	assert(if_def_UI32[id]);
	return states_UI32[id];
}

void* States::getptr (t_PTR id) {
	assert(id < PTR_LASTID);
	assert(if_def_PTR[id]);
	return states_PTR[id];
}

void States::setdouble(t_DOUBLE id, double d) {
	assert(id < DOUBLE_LASTID);
	if_def_DOUBLE[id] = true;
	states_DOUBLE[id] = d;
}

void States::setuint16(t_UI16 id, uint16 d) {
	assert(id < UI16_LASTID);
	if_def_UI16[id] = true;
	states_UI16[id] = d;
}

void States::setuint32(t_UI32 id, uint32 d) {
	assert(id < UI32_LASTID);
	if_def_UI32[id] = true;
	states_UI32[id] = d;
}

void States::setptr(t_PTR id, void* d) {
	assert(id < PTR_LASTID);
	if_def_PTR[id] = true;
	states_PTR[id] = d;
}

bool States::ifdefdouble (t_DOUBLE id) {
		assert(id < DOUBLE_LASTID);
		return if_def_DOUBLE[id];
}

bool States::ifdefuint16 (t_UI16 id) {
		assert(id < UI16_LASTID);
		return if_def_UI16[id];
}

bool States::ifdefuint32 (t_UI32 id) {
		assert(id < UI32_LASTID);
		return if_def_UI32[id];
}

bool States::ifdefptr (t_PTR id) {
		assert(id < PTR_LASTID);
		return if_def_PTR[id];
}

void States::undefinedouble (t_DOUBLE id) {
		assert(id < DOUBLE_LASTID);
		if_def_DOUBLE[id] = false;
}
void States::undefineuint16 (t_UI16 id) {
		assert(id < UI16_LASTID);
		if_def_UI16[id] = false;
}

void States::undefineuint32 (t_UI32 id) {
		assert(id < UI32_LASTID);
		if_def_UI32[id] = false;
}

void States::undefineptr (t_PTR id) {
		assert(id < PTR_LASTID);
		if_def_PTR[id] = false;
}

} // namespace nla3d

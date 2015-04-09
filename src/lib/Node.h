// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include "math/Vec.h"
#include "Dof.h"

namespace nla3d {

//class Node represents spatial 3D node
class Node {
public:
	Node();
	~Node() { }
  // working with Dofs
	static uint16 getNumberOfDofs();	
  static void registerDofType(Dof::dofType type);
	static Dof::dofType getDofType (uint16 dofIndex);
  static uint16 getDofIndex (Dof::dofType dof);
  static bool isDofUsed (Dof::dofType dof);

	//in-out operation: (rudiment actually..)
	void display (uint32 nn);
	std::string toString();

  math::Vec<3> pos;

	friend class Element;
private:
  // vector for mapping dofType into number of Dof in current task
  static std::vector<uint16> dofNumberList;
	static uint16 numberOfDofs;
};

inline uint16 Node::getNumberOfDofs() {
  return numberOfDofs;		
}

inline uint16 Node::getDofIndex (Dof::dofType dof) {
  return dofNumberList[dof];
}


inline bool Node::isDofUsed (Dof::dofType dof) {
  if (dofNumberList[dof] != Dof::UNDEFINED) {
    return true;
  } else {
    return false;
  }
}

inline void Node::display (uint32 nn) {
  LOG(INFO) << "N " << nn << " : " << pos;
}

inline std::string Node::toString()	{
  std::string str;
  char buff[100];
  sprintf_s(buff,100,"%f %f %f",pos[0], pos[1], pos[2]);
  str+=buff;
  return str; //TODO: do it easy
}

} // namespace nla3d 

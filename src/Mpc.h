// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include <list>
#include "sys.h"
#include "Dof.h"

namespace nla3d {

class FEStorage;

class MpcTerm {
public:
	MpcTerm() : node(0), node_dof(Dof::UNDEFINED), coef(1.0) {

  }

	MpcTerm(int32 n, Dof::dofType dof, double coef = 1.0) : node(n), node_dof(dof), coef(coef) {
  
  }

	int32 node;
  Dof::dofType node_dof;
	double coef;
};


class Mpc {
  public:
  bool isUpdatable;
  bool isSinglePointConstraint;
  //virtual void update(double time, double dtime) = 0;
  void print (std::ostream& out);
	std::list<MpcTerm> eq;
	double b;
  //virtual void pre() { };
};


class MpcCollection {
  public:
  virtual void update() = 0;
  virtual void pre() = 0;
  void printEquations (std::ostream& out);
  void registerMpcsInStorage();
  std::vector<Mpc*> collection;
  FEStorage* storage;
};


class RigidBodyMpc : public MpcCollection {
  public:
  RigidBodyMpc();
  ~RigidBodyMpc();
  uint32 masterNode;
  
  void pre ();
  void update ();
  std::vector<uint32> slaveNodes;
  std::vector<Dof::dofType> dofs;
};

// old fashioned class BC and its derivatives
class BC {
public:
	BC () : isUpdatetable(false) {  }
	//virtual void apply(FEStorage *storage)=0;
  // does this BC need to be updated every step (new eq. every step)
	bool isUpdatetable; //нужно ли обновлять на каждом шаге решения
};

class BC_dof_constraint : public BC {
public:
	BC_dof_constraint () : node(0), node_dof(Dof::UNDEFINED), value(0.0)
	{	}
	int32 node;
  Dof::dofType node_dof;
	double value;
	//void apply(FEStorage *storage);
};

class BC_dof_force : public BC {
public:
	BC_dof_force () : node(0), node_dof(Dof::UNDEFINED), value(0.0)
	{	}
	int32 node;
  Dof::dofType node_dof;
	double value;
	//void apply(FEStorage *storage);
};


//class MpcFiniteRotation : public Mpc {
//  //unit axis of rotation
//  math::Vec<3> n;
//  //zero point of rotation 
//  math::Vec<3> o; 
//  // speed of rotation
//  double theta;
//
//  uint32 masterNode;
//  
//  vector<uint32> slaveNodes;
//  vector<Dof::dofType> dofs;
//};


// MPC equation like
// coef1 * dof1 + coef2 * dof2 + ... + coefn * dofn - b = 0
//class BC_MPC : public BC {
//public:
//	std::list<MPC_token> eq;
//	double b;
//	// void apply(FEStorage *storage);
//  // virtual void update(const double time) { };
//};


} // namespace nal3d

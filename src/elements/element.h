// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include <vector>
#include <math.h>
#include "query.h"
#include "Node.h"
#include "math\Vec.h"
#include "math\Mat.h"

namespace nla3d {

//pre-defines
class Material;
class FEStorage;


//----------------------------------------------//
//--------------------- Face -------------------//
//----------------------------------------------//
class Face
{
public:
	friend class Element;
	Face() { }
	Face(uint32 n1, uint32 n2, uint32 n3, uint32 n4=0)
	{
		nodes[0] = n1;
		nodes[1] = n2;
		nodes[2] = n3;
		nodes[3] = n4;
	}
	uint32 nodes[4]; //only QUAD and TRIANGLE are supported
	friend std::ostream& operator<< (std::ostream& stream,const Face& obj);
private:
	static uint16 num_nodes_per_face; //only class Element can change it
};


//----------------------------------------------//
//-------------------- Element------------------//
//----------------------------------------------//
//Element abstract class
class Element {
public: 
	Element () : nodes(NULL) {
    // initialise static variables (just once)
    if (dofNumberList.size() == 0) {
      dofNumberList.assign(Dof::numberOfDofTypes, Dof::UNDEFINED);
      numberOfDofs = 0;
    }
	}
	~Element() {
		if (nodes) delete[] nodes;
    nodes = NULL;
	}
	static uint16 n_nodes();
	static uint16 n_dim();
	static uint16 n_dofs();
	uint32& getNodeNumber (uint16 num);
	static uint16 n_int();
	static void set_n_int(uint16 _nint);// нельзя вызывать после выполнения функции pre() (начало решения)
	static uint16 get_central_gp();
	static uint16 n_face();
	static uint16 n_nodes_per_face();
	static Face& get_face(uint16 fn);
  static FEStorage& getStorage();

  // heart of the element class
	virtual void pre()=0;
	virtual void build()=0;
	virtual void update()=0;
	virtual void getScalar(double& scalar, query::scalarQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0)=0;
	virtual void getVector(double* vector, query::vectorQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0)=0;
	virtual void getTensor(math::MatSym<3>& tensor, query::tensorQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0)=0;

	Element& operator= (const Element& from);

  uint32 getElNum();

	//in-out operation:
	void read_from_stream (std::istream& str);
	void print (std::ostream& out);


  // working with Dofs
	static uint16 getNumberOfDofs();	
  static void registerDofType(Dof::dofType type);
	static Dof::dofType getDofType (uint16 dofIndex);
  static uint16 getDofIndex (Dof::dofType dof);
  static bool isDofUsed (Dof::dofType dof);

  friend class FEStorage;
protected:

	void change_face_nodes_num (uint16 num);

	static uint16 number_of_nodes;
	static uint16 number_of_dimensions;
	static std::vector<Face> faces;
	static uint16 number_of_integration_points;	// number of int points per coordinate
	uint32 *nodes;
  uint32 elNum;
  static FEStorage* storage;


  // vector for mapping dofType into number of Dof in current task
  static std::vector<uint16> dofNumberList;
	static uint16 numberOfDofs;
};

//----------------------------------------------//
//-----------------Element_SOLID8---------------//
//----------------------------------------------//
//Element geometry class (on this step we describe how many nodes and node's DOFs)
class Element_SOLID8 : public Element {
public:
	Element_SOLID8 ()
	{
		Element::number_of_nodes = 8;
		number_of_dimensions = 3;
		nodes = new uint32[Element::number_of_nodes];
		Element::change_face_nodes_num(4);
		if (faces.empty())
		{
			faces.push_back(Face(1,2,3,4));
			faces.push_back(Face(8,7,6,5));
			faces.push_back(Face(7,3,2,6));
			faces.push_back(Face(8,5,1,4));
			faces.push_back(Face(4,3,7,8));
			faces.push_back(Face(1,5,6,2));
		}
	}
	Element_SOLID8& operator= (const Element_SOLID8& from)
	{
		Element::operator= (from);
		return *this;
	}
};

//----------------------------------------------//
//-----------------Element_PLANE4---------------//
//----------------------------------------------//
class Element_PLANE4 : public Element
{
public:
	Element_PLANE4 ()
	{
		//change_node_dofs_num(2);
		Element::number_of_nodes = 4;
		number_of_dimensions = 2;
		nodes = new uint32[Element::number_of_nodes];
		Element::change_face_nodes_num(4);
		if (faces.empty())
		{
			faces.push_back(Face(1,2,3,4));
		}
	}
	Element_PLANE4& operator= (const Element_PLANE4& from)
	{
		Element::operator= (from);
		return *this;
	}
};



inline uint16 Element::n_nodes() {
	return number_of_nodes;
}
inline uint16 Element::n_dim() {
	return number_of_dimensions;
}
inline uint16 Element::n_dofs() {
	return numberOfDofs;
}
inline uint32& Element::getNodeNumber (uint16 num) { //TODO: нужна ли ссылка тут?
	assert(num < n_nodes());
	assert(nodes);
	return nodes[num];
}
inline uint16 Element::n_int() {
	return number_of_integration_points;
}

inline void Element::set_n_int(uint16 _nint) { //TODO: make more safety: разрешить менять до операции pre()
	number_of_integration_points = _nint; 
}

inline uint16 Element::get_central_gp() {
	return (uint16) npow(n_int(), n_dim())/2; //TODO: is it appropriate?
}

inline uint16 Element::n_face() {
	return static_cast<uint16> (faces.size());
}

inline uint16 Element::n_nodes_per_face() {
	return Face::num_nodes_per_face;
}

inline Face& Element::get_face(uint16 fn) {
	assert(fn < n_face());
	return faces[fn];
}

inline void Element::change_face_nodes_num (uint16 num) {
	Face::num_nodes_per_face = num;
}

inline bool Element::isDofUsed (Dof::dofType dof) {
  if (dofNumberList[dof] != Dof::UNDEFINED) {
    return true;
  } else {
    return false;
  }
}

inline FEStorage& Element::getStorage() {
  return *storage;
}

inline uint32 Element::getElNum() {
  return elNum;
}

inline uint16 Element::getNumberOfDofs() {
  return numberOfDofs;		
}

inline uint16 Element::getDofIndex (Dof::dofType dof) {
  return dofNumberList[dof];
}

} // namespace nla3d

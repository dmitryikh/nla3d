#pragma once
#include "sys.h"
#include <vector>
#include <math.h>
#include "math\Vec.h"
#include "math\Mat.h"
#include "FE_Storage.h"

//pre-defines
class Material;
class FE_Storage_Interface;

enum Dof_Type
{
	UX,
	UY,
	UZ,
	ROTX,
	ROTY,
	ROTZ,
	HYDRO_PRESSURE,
	UNKNOWN
};

// class for Degree of Freedom informations 
// (is it predefined by BC, its number in equation system, ..)
class Dof
{
public:
	Dof() : eq_number(0), is_constrained(false)
	{
	}
	uint32 eq_number; // 0 - not use
	bool is_constrained; // is this DOF defined by B.C.
};

//class Node represents spatial 3D node
class Node
{
public:
	Node()
	{	}
	~Node()
	{	}
	static uint16 n_dofs()
	{
		return number_of_dofs;		
	}
	
	static Dof_Type dof_type(uint16 dof)
	{
		assert(dof < dof_types.size());
		return dof_types[dof];
	}
	//in-out operation:
	void display (uint32 nn)
	{
		echo("N %d: %f\t%f\t%f", nn, pos[0], pos[1], pos[2]);
	}
	string toString()
	{
		string str;
		char buff[100];
		sprintf_s(buff,100,"%f %f %f",pos[0], pos[1], pos[2]);
		str+=buff;
		return str; // do it easy
	}
	void read_from_stream(istream &str)
	{
		str >> pos[0];
		str >> pos[1];
		str >> pos[2]; //TODO: process wrong input
	}

	Vec<3> pos;

	friend class Element;
private:
	static uint16 number_of_dofs; //only class Element can change it
	static vector<Dof_Type> dof_types; //only class Element can change it
};

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

//Element abstract class
class Element {
public: 
  enum elTypes {
    NOT_DEFINED,
    PLANE41,
    SOLID81,
    LAST
  };

  static const char* const elTypeLabels[]={"UNDEFINED",
    "PLANE41",
    "SOLID81",
  };

	Element () : nodes(NULL) {
	}
	~Element() {
		if (nodes) delete[] nodes;
    nodes = NULL;
	}
	static uint16 n_nodes();
	static uint16 n_dim();
	static uint16 n_dofs();
	uint32& node_num (uint16 num);
	static uint16 n_int();
	static void set_n_int(uint16 _nint);// нельзя вызывать после выполнения функции pre() (начало решения)
	static uint16 get_central_gp();
	static uint16 n_face();
	static uint16 n_nodes_per_face();
	static Face& get_face(uint16 fn);
	static Dof_Type dof_type(uint16 dof);
  static FE_Storage_Interface& getStorage();

	template <uint16 el_dofs_num>
	void assemble (const Mat<el_dofs_num,el_dofs_num> &Ke, const Vec<el_dofs_num> &Qe);
	template <uint16 dimM, uint16 dimN>
	void assemble2(MatSym<dimM> &Kuu, Mat2<dimM,dimM> &Kup, Mat2<dimN,dimN> &Kpp, Vec<dimM> &Fu, Vec<dimN> &Fp);
	template <uint16 dimM>
	void assemble3(MatSym<dimM> &Kuu, Vec<dimM> &Kup, double Kpp, Vec<dimM> &Fu, double Fp);

  // hear of the element class
	virtual void pre()=0;
	virtual void build()=0;
	virtual void update()=0;
	virtual void getScalar(double& scalar, el_component code, uint16 gp = GP_MEAN, const double scale = 1.0)=0;
	virtual void getTensor(MatSym<3>& tensor, el_tensor code, uint16 gp = GP_MEAN, const double scale = 1.0)=0;

	Element& operator= (const Element& from);

  uint32 getElNum();

  static elTypes elName2elType (string elName); 
  static void createElements (elTypes elId, const uint32 n, void** ptr); 

	//in-out operation:
	void read_from_stream (istream &str);
	void display (uint32 en);
	string toString();

  friend class FE_Storage_Interface;
protected:

	void change_node_dofs_num (uint16 ndof,...);
	void change_el_dofs_num (uint16 ndof, ...);
	void change_face_nodes_num (uint16 num);

	static uint16 number_of_nodes;
	static uint16 number_of_dimensions;
	static vector<Face> faces;
	static uint16 number_of_integration_points;	// number of int points per coordinate
	static uint16 number_of_dofs;
	static vector<Dof_Type> dof_types;
	uint32 *nodes;
  uint32 elNum;
  static FE_Storage_Interface* storage;
};


//Element geometry class (on this step we describe how many nodes and node's DOFs)
class Element_8N_3D : public Element
{
public:
	Element_8N_3D ()
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
	Element_8N_3D& operator= (const Element_8N_3D& from)
	{
		Element::operator= (from);
		return *this;
	}
};

class Element_4N_2D : public Element
{
public:
	Element_4N_2D ()
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
	Element_4N_2D& operator= (const Element_4N_2D& from)
	{
		Element::operator= (from);
		return *this;
	}
};



inline uint16 Element::n_nodes()
{
	return number_of_nodes;
}
inline uint16 Element::n_dim()
{
	return number_of_dimensions;
}
inline uint16 Element::n_dofs()
{
	return number_of_dofs;
}
inline uint32& Element::node_num (uint16 num) //TODO: нужна ли ссылка тут?
{
	assert(num < n_nodes());
	assert(nodes);
	return nodes[num];
}
inline uint16 Element::n_int()
{
	return number_of_integration_points;
}
inline void Element::set_n_int(uint16 _nint) //TODO: make more safety: разрешить менять до операции pre()
{
	number_of_integration_points = _nint; 
}
inline uint16 Element::get_central_gp()
{
	return (uint16) npow(n_int(), n_dim())/2; //TODO: is it appropriate?
}
inline uint16 Element::n_face()
{
	return faces.size();
}
inline uint16 Element::n_nodes_per_face()
{
	return Face::num_nodes_per_face;
}
inline Face& Element::get_face(uint16 fn)
{
	assert(fn < n_face());
	return faces[fn];
}
inline void Element::change_face_nodes_num (uint16 num)
{
	Face::num_nodes_per_face = num;
}
inline Dof_Type Element::dof_type(uint16 dof)
{
	assert(dof < dof_types.size());
	return dof_types[dof];
}

inline FE_Storage_Interface& Element::getStorage() {
  return *storage;
}

inline uint32 Element::getElNum() {
  return elNum;
}


template <uint16 el_dofs_num>
void Element::assemble (const Mat<el_dofs_num,el_dofs_num> &Ke, const Vec<el_dofs_num> &Qe)
{
	uint16 eds = Element::n_nodes()*Node::n_dofs(); // el's dofs start number
	assert(el_dofs_num == eds+Element::n_dofs());
	// upper diagonal process for nodal dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 j=i; j < Element::n_nodes(); j++)
			for (uint16 di=0; di < Node::n_dofs(); di++)
				for (uint16 dj=0; dj < Node::n_dofs(); dj++)
					if ((i==j) && (dj>di)) continue;
					else
						storage->Kij_add(nodes[i],di,nodes[j],dj, Ke[i*Node::n_dofs()+di][j*Node::n_dofs()+dj]);
	//upper diagonal process for nodes-el dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for(uint16 di=0; di < Node::n_dofs(); di++)
			for (uint16 dj=0; dj < Element::n_dofs(); dj++)
				storage->Kij_add(nodes[i],di, -(int32)getElNum(), dj, Ke[i*Node::n_dofs()+di][eds+dj]);
	//upper diagonal process for el-el dofs
	for (uint16 di=0; di < Element::n_dofs(); di++)
		for (uint16 dj=di; dj < Element::n_dofs(); dj++)
			storage->Kij_add(-(int32)getElNum(), di, -(int32)getElNum(), dj,  Ke[eds+di][eds+dj]);

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
			storage->Fi_add(nodes[i],di, Qe[i*Node::n_dofs()+di]);

	for (uint16 di=0; di < Element::n_dofs(); di++)
		storage->Fi_add(-(int32)getElNum(), di, Qe[eds+di]);
}

template <uint16 dimM, uint16 dimN>
void Element::assemble2(MatSym<dimM> &Kuu, Mat2<dimM,dimM> &Kup, Mat2<dimN,dimN> &Kpp, Vec<dimM> &Fu, Vec<dimN> &Fp) 
{
	assert (Element::n_nodes()*Node::n_dofs() == dimM);
	assert (Element::n_dofs() == dimN);
	double *Kuu_p = Kuu.ptr();
	double *Kup_p = Kup.ptr();
	double *Kpp_p = Kpp.ptr();
	double *Fu_p = Fu_p.ptr();
	double *Fp_p = Fp.ptr();                               

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
			for (uint16 j=i; j < Element::n_nodes(); j++)
				for (uint16 dj=di; dj < Node::n_dofs(); dj++)
				{
						storage->Kij_add(nodes[i],di,nodes[j],dj, *Kuu_p);
						Kuu_p++;
				}
	//upper diagonal process for nodes-el dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for(uint16 di=0; di < Node::n_dofs(); di++)
			for (uint16 dj=0; dj < Element::n_dofs(); dj++)
			{
				storage->Kij_add(nodes[i],di, -(int32)getElNum(), dj, *Kup_p);
				Kup_p++;
			}
	//upper diagonal process for el-el dofs
	for (uint16 di=0; di < Element::n_dofs(); di++)
		for (uint16 dj=di; dj < Element::n_dofs(); dj++)
		{
			storage->Kij_add(-(int32)getElNum(), di, -(int32)getElNum(), dj,  *Kpp_p);
			Kpp_p++;
		}

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
		{
			storage->Fi_add(nodes[i],di, *Fu_p);
			Fu_p++;
		}

	for (uint16 di=0; di < Element::n_dofs(); di++)
	{
		storage->Fi_add(-(int32)getElNum(), di, *Fp_p);
		Fp_p++;
	}
}

template <uint16 dimM>
void Element::assemble3(MatSym<dimM> &Kuu, Vec<dimM> &Kup, double Kpp, Vec<dimM> &Fu, double Fp) 
{
	assert (Element::n_nodes()*Node::n_dofs() == dimM);
	assert (Element::n_dofs() == 1);
	double *Kuu_p = Kuu.ptr();
	double *Kup_p = Kup.ptr();
	double *Fu_p = Fu.ptr();

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
		for (uint16 j=i; j < Element::n_nodes(); j++)
			
				for (uint16 dj=0; dj < Node::n_dofs(); dj++)
				{
					if ((i==j) && (dj<di)) continue;
					else
					{
						storage->Kij_add(nodes[i],di,nodes[j],dj, *Kuu_p);
						Kuu_p++;
					}
				}
	//upper diagonal process for nodes-el dofs
	for (uint16 i=0; i < Element::n_nodes(); i++)
		for(uint16 di=0; di < Node::n_dofs(); di++)
			for (uint16 dj=0; dj < Element::n_dofs(); dj++)
			{
				storage->Kij_add(nodes[i],di, -(int32)getElNum(), dj, *Kup_p);
				Kup_p++;
			}
	//upper diagonal process for el-el dofs
			storage->Kij_add(-(int32)getElNum(), 0, -(int32)getElNum(), 0,  Kpp);

	for (uint16 i=0; i < Element::n_nodes(); i++)
		for (uint16 di=0; di < Node::n_dofs(); di++)
		{
			storage->Fi_add(nodes[i],di, *Fu_p);
			Fu_p++;
		}
		storage->Fi_add(-(int32)getElNum(), 0, Fp);
}

  

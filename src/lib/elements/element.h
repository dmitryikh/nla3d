// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include <vector>
#include <math.h>
#include <initializer_list>
#include <Eigen/Dense>

#include "query.h"
#include "Node.h"
#include "FEStorage.h"
#include "math/Vec.h"
#include "math/Mat.h"

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
  Element () : nodes(nullptr) { }
  virtual ~Element() {
    if (nodes) {
      delete[] nodes;
      nodes = nullptr;
    }
  }
  static uint16 n_nodes();
  static uint16 n_dim();
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
  virtual void getScalar(double& scalar, query::scalarQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0);
  virtual void getVector(double* vector, query::vectorQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0);
  virtual void getTensor(math::MatSym<3>& tensor, query::tensorQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0);

  Element& operator= (const Element& from);

  uint32 getElNum();

  //in-out operation:
  void print (std::ostream& out);

  // test implementation
  template <uint16 dimM>
  void assemble (math::MatSym<dimM> &Ke, std::initializer_list<Dof::dofType> _nodeDofs);
  void assemble (Eigen::Ref<Eigen::MatrixXd> Ke, std::initializer_list<Dof::dofType> _nodeDofs);

  friend class FEStorage;
protected:

  void change_face_nodes_num (uint16 num);

  static uint16 number_of_nodes;
  static uint16 number_of_dimensions;
  static std::vector<Face> faces;
  static uint16 number_of_integration_points; // number of int points per coordinate
  uint32 *nodes;
  uint32 elNum;
  static FEStorage* storage;
};

//----------------------------------------------//
//-----------------Element_SOLID8---------------//
//----------------------------------------------//
//Element geometry class (on this step we describe how many nodes and node's DOFs)
class Element_SOLID8 : public Element {
public:
  Element_SOLID8 () {
    Element::number_of_nodes = 8;
    number_of_dimensions = 3;
    nodes = new uint32[Element::number_of_nodes];
    Element::change_face_nodes_num(4);
    if (faces.empty()) {
      faces.push_back(Face(1,2,3,4));
      faces.push_back(Face(8,7,6,5));
      faces.push_back(Face(7,3,2,6));
      faces.push_back(Face(8,5,1,4));
      faces.push_back(Face(4,3,7,8));
      faces.push_back(Face(1,5,6,2));
    }
  }
  Element_SOLID8& operator= (const Element_SOLID8& from) {
    Element::operator= (from);
    return *this;
  }
};

//----------------------------------------------//
//-----------------Element_PLANE4---------------//
//----------------------------------------------//
class Element_PLANE4 : public Element {
public:
  Element_PLANE4 () {
    //change_node_dofs_num(2);
    Element::number_of_nodes = 4;
    number_of_dimensions = 2;
    nodes = new uint32[Element::number_of_nodes];
    Element::change_face_nodes_num(4);
    if (faces.empty()) {
      faces.push_back(Face(1,2,3,4));
    }
  }
  Element_PLANE4& operator= (const Element_PLANE4& from) {
    Element::operator= (from);
    return *this;
  }
};

template <uint16 dimM>
void Element::assemble (math::MatSym<dimM> &Ke, std::initializer_list<Dof::dofType> _nodeDofs) {
  assert (nodes != NULL);
  double* Ke_p = Ke.ptr();
  std::vector<Dof::dofType> nodeDof(_nodeDofs);
  uint16 dim = static_cast<uint16> (_nodeDofs.size());  // TODO: not true!
  assert (Element::n_nodes() * dim == dimM);

  for (uint16 i=0; i < Element::n_nodes(); i++) {
    for (uint16 di=0; di < dim; di++) {
      for (uint16 j=i; j < Element::n_nodes(); j++) {
        for (uint16 dj=0; dj < dim; dj++) {
          if ((i==j) && (dj<di)) {
            continue;
          } else {
            storage->Kij_add(nodes[i], nodeDof[di], nodes[j], nodeDof[dj], *Ke_p);
            Ke_p++;
          }
        }
      }
    }
  }
}

// may be useful way to build a global stiffnes matrix from element's stiff. matrix
//double* Ke_p = Ke.ptr();
//Dof::dofType dofi;
//uint32 nodei;
//Dof::dofType dofj;
//uint32 nodej;
//for (uint16 i = 0; i < Element::n_dofs(); i++) {
//  getDofByRowNumber(i, &nodei, &dofi); 
//  // this is for diagonal elements (we economy with getDofByRowNumber calls here)
//  storage->Kij_add(nodei, dofi, nodei, dofi, *Ke_p);
//  Ke_p++;
//  for (uint16 j = i + 1; j < Element::n_dofs(); j++) {
//    getDofByRowNumber(j, &nodej, &dofj); 
//    storage->Kij_add(nodei, dofi, nodej, dofj, *Ke_p);
//    Ke_p++;
//  }
//}

inline uint16 Element::n_nodes() {
  return number_of_nodes;
}


inline uint16 Element::n_dim() {
  return number_of_dimensions;
}
// & is used here because this function is called such this:
// el->getNodeNumber(0) = 1234;
inline uint32& Element::getNodeNumber (uint16 num) {
  assert(num < n_nodes());
  assert(nodes);
  return nodes[num];
}

inline uint16 Element::n_int() {
  return number_of_integration_points;
}

// TODO: how to avoid changing of number_of_integration_points after pre() command?
inline void Element::set_n_int(uint16 _nint) {
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


inline FEStorage& Element::getStorage() {
  return *storage;
}

inline uint32 Element::getElNum() {
  return elNum;
}

} // namespace nla3d

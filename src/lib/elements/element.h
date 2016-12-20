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


// Element shapes are taken from vtk-file-formats.pdf
enum ElementShape {
  VERTEX = 0,
  LINE,
  TRIANGLE,
  QUAD,
  TETRA,
  HEXAHEDRON,
  WEDGE,
  PYRAMID,
  QUADRARIC_EDGE,
  QUADRATIC_TRIANGLE,
  QUADRATIC_QUAD,
  QUADRATIC_TETRA,
  QUADRARIC_HEXAHEDRON,
  UNDEFINED
};

// number of dimensions in shape
static const uint16 _shape_dim[] = {
  0,  // VERTEX
  1,  // LINE
  2,  // TRIANGLE
  2,  // QUAD
  3,  // TETRA
  3,  // HEXAHEDRON
  3,  // WEDGE
  3,  // PYRAMID
  2,  // QUADRARIC_EDGE
  2,  // QUADRATIC_TRIANGLE
  2,  // QUADRATIC_QUAD
  3,  // QUADRATIC_TETRA
  3  // QUADRARIC_HEXAHEDRON
};


// number of nodes in shape
static const uint16 _shape_nnodes[] = {
  1,  // VERTEX
  2,  // LINE
  3,  // TRIANGLE
  4,  // QUAD
  4,  // TETRA
  8,  // HEXAHEDRON
  6,  // WEDGE
  5,  // PYRAMID
  3,  // QUADRARIC_EDGE
  6,  // QUADRATIC_TRIANGLE
  8,  // QUADRATIC_QUAD
  10, // QUADRATIC_TETRA
  20  // QUADRARIC_HEXAHEDRON
};


//Element abstract class
class Element {
  public: 
    Element ();
    virtual ~Element();

    uint32 getElNum();
    uint16 getNNodes();
    uint16 getDim();
    ElementShape getShape();
    ElementType getType();
    uint32& getNodeNumber (uint16 num);
    FEStorage& getStorage();
    uint16 getIntegrationPoints();
    void setIntegrationPoints(uint16 _nint); // нельзя вызывать после выполнения функции pre() (начало решения)

    // heart of the element class
    virtual void pre()=0;
    virtual void build()=0;
    virtual void update()=0;
    virtual void getScalar(double& scalar, query::scalarQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0);
    virtual void getVector(double* vector, query::vectorQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0);
    virtual void getTensor(math::MatSym<3>& tensor, query::tensorQuery code, uint16 gp = query::GP_MEAN, const double scale = 1.0);

    Element& operator= (const Element& from);

    //in-out operation:
    void print (std::ostream& out);

    // may be useful for derived particular class
    template <uint16 dimM>
    void assemble (math::MatSym<dimM> &Ke, std::initializer_list<Dof::dofType> _nodeDofs);
    void assemble (Eigen::Ref<Eigen::MatrixXd> Ke, std::initializer_list<Dof::dofType> _nodeDofs);

    friend class FEStorage;
  protected:

    ElementType type = ElementType::UNDEFINED;
    ElementShape shape = ElementShape::UNDEFINED;
    uint16 nOfIntPoints = 0; // number of int points overall
    uint32 elNum = 0;
    uint32 *nodes = nullptr;
    FEStorage* storage = nullptr;
};


//Element geometry class
class ElementLINE : public Element {
  public:
    ElementLINE() {
      shape = ElementShape::LINE;
      nodes = new uint32[getNNodes()];
    }

    ElementLINE& operator= (const ElementLINE& from) {
      Element::operator= (from);
      return *this;
    }
};


class ElementQUAD : public Element {
  public:
    ElementQUAD () {
      shape = ElementShape::QUAD;
      nodes = new uint32[getNNodes()];
    }

    ElementQUAD& operator= (const ElementQUAD& from) {
      Element::operator= (from);
      return *this;
    }
};


class ElementTETRA : public Element {
  public:
    ElementTETRA() {
      shape = ElementShape::TETRA;
      nodes = new uint32[getNNodes()];
    }

    ElementTETRA& operator= (const ElementTETRA& from) {
      Element::operator= (from);
      return *this;
    }
};


class ElementHEXAHEDRON : public Element {
  public:
    ElementHEXAHEDRON() {
      shape = ElementShape::HEXAHEDRON;
      nodes = new uint32[getNNodes()];
    }

    ElementHEXAHEDRON& operator= (const ElementHEXAHEDRON& from) {
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
  assert (getNNodes() * dim == dimM);

  for (uint16 i=0; i < getNNodes(); i++) {
    for (uint16 di=0; di < dim; di++) {
      for (uint16 j=i; j < getNNodes(); j++) {
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


inline uint16 Element::getNNodes() {
  return _shape_nnodes[shape];
}


inline uint16 Element::getDim() {
  return _shape_dim[shape];
}


inline ElementShape Element::getShape() {
  return shape;
}


// & is used here because this function is called such this:
// el->getNodeNumber(0) = 1234;
inline uint32& Element::getNodeNumber (uint16 num) {
  assert(num < getNNodes());
  assert(nodes);
  return nodes[num];
}


inline uint16 Element::getIntegrationPoints() {
  return nOfIntPoints;
}


// TODO: how to avoid changing of nOfIntPoints after pre() command?
inline void Element::setIntegrationPoints(uint16 _nint) {
  nOfIntPoints = _nint; 
}


inline FEStorage& Element::getStorage() {
  return *storage;
}


inline uint32 Element::getElNum() {
  return elNum;
}

} // namespace nla3d

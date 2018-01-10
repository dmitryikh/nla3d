%module nla3d
%{
#define SWIG_FILE_WITH_INIT
#include "math/Vec.h"
#include "math/Mat.h"
#include "FEStorage.h"
#include "PostProcessor.h"
#include "VtkProcessor.h"
#include "FESolver.h"
#include "FEReaders.h"
#include "elements/query.h"
#include "elements/TRUSS3.h"
#include "elements/TETRA0.h"
#include "elements/QUADTH.h"
using namespace nla3d;
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"

%rename($ignore, %$isfunction) ""; // ignore all functions
%rename($ignore, %$isclass) ""; // ignore all classes
%ignore "UNDEFINED";
%include "sys.h"
%rename("%s", %$isfunction) "";
%rename("%s", %$isclass) "";

// global rules for c++ exceptions handling
%exception {
  try {
    $action
  }
  catch(std::range_error& e) {
      PyErr_SetString(PyExc_IndexError, e.what());
      SWIG_fail;
  }
  catch(std::bad_cast& e) {
      PyErr_SetString(PyExc_IndexError, e.what());
      SWIG_fail;
  }
}

// Vec.h
%ignore nla3d::math::dVec::dVec(nla3d::math::dVec &&);
%ignore nla3d::math::dVec::operator[](uint32 _n);
%ignore nla3d::math::dVec::operator[](uint32 _n) const;
%ignore nla3d::math::dVec::operator=;
%ignore nla3d::math::Vec::operator[](uint16 n);
%ignore nla3d::math::Vec::operator[](uint16 n) const;
%ignore nla3d::math::Vec::operator=;

%include "math/Vec.h"

// Vec<dimM>
%define VEC_EXTEND(N)
%template(Vec ## N) nla3d::math::Vec<N>;
%extend nla3d::math::Vec<N> {
  %pythoncode { 
    def __str__(self): 
      _list = []
      for i in range(self.size()):
        _list.append(self[i])
      return str(_list)
  } 
  double __getitem__(uint16 _n ) {
    if (_n >= $self->size())
      throw std::range_error("Index out of range");
    return $self->operator[](_n);
  }
  void __setitem__(uint16 _n, double val) {
    if (_n >= $self->size())
      throw std::range_error("Index out of range");
    $self->operator[](_n) = val;
  }
}
%enddef

VEC_EXTEND(2)
VEC_EXTEND(3)
VEC_EXTEND(4)
VEC_EXTEND(5)
VEC_EXTEND(6)
VEC_EXTEND(7)

//dVec
%extend nla3d::math::dVec {
  %pythoncode { 
    def __str__(self): 
      _list = []
      for i in range(self.size()):
        _list.append(self[i])
      return str(_list)
  } 
  double __getitem__(uint32 _n ) {
    if (_n >= $self->size())
      throw std::range_error("Index out of range");
    return $self->operator[](_n);
  }
  void __setitem__(uint32 _n, double val) {
    if (_n >= $self->size())
      throw std::range_error("Index out of range");
    $self->operator[](_n) = val;
  }
}

// Mat.h
%ignore nla3d::math::dMat_interface;
%ignore nla3d::math::Mat::cross_cut;
%include "math/Mat.h"

%define MAT_EXTEND(M, N)
%template(Mat ## M ## N) nla3d::math::Mat<M, N>;
%extend nla3d::math::Mat<M, N> {
  %pythoncode { 
    def __str__(self):
      return str(self.to_list())

    def __setitem__(self, pos, val):
      i, j = pos
      self._setitem(i, j, val)

    def __getitem__(self, pos):
      i, j = pos
      return self._getitem(i, j)

    def from_list(self, list):
      assert(len(list) == self.dM())
      for i, ll in enumerate(list):
        assert(len(ll) == self.dN())
        for j, val in enumerate(ll):
          self[i, j] = val

    def to_list(self):
      list = []
      for i in range(self.dM()):
        ll = []
        for j in range(self.dN()):
          ll.append(self[i, j])
        list.append(ll)
      return list
  } 

  double _getitem(uint16 i, uint16 j) {
    if (i >= $self->dM() || j >= $self->dN())
      throw std::range_error("Index out of range");
    return $self->operator[](i)[j];
  }

  void _setitem(uint16 i, uint16 j, double val) {
    if (i >= $self->dM() || j >= $self->dN())
      throw std::range_error("Index out of range");
    $self->operator[](i)[j] = val;
  }

  Mat<M, N> copy() {
    return Mat<M, N>(*$self);
  }
}
%enddef

MAT_EXTEND(3,3)

%extend nla3d::math::dMat {
  %pythoncode {
    def __str__(self):
      return str(self.to_list())

    def __setitem__(self, pos, val):
      i, j = pos
      self._setitem(i, j, val)

    def __getitem__(self, pos):
      i, j = pos
      return self._getitem(i, j)

    @staticmethod
    def from_list(list):
      assert(len(list) > 0)
      assert(len(list[0]) > 0)
      M, N = len(list), len(list[0])
      mat = dMat(M, N)
      for i, ll in enumerate(list):
        assert(len(ll) == mat.dN())
        for j, val in enumerate(ll):
          mat[i, j] = val
      return mat

    def to_list(self):
      list = []
      for i in range(self.dM()):
        ll = []
        for j in range(self.dN()):
          ll.append(self[i, j])
        list.append(ll)
      return list
  }

  double _getitem(uint16 i, uint16 j) {
    if (i >= $self->dM() || j >= $self->dN())
      throw std::range_error("Index out of range");
    return $self->operator[](i)[j];
  }

  void _setitem(uint16 i, uint16 j, double val) {
    if (i >= $self->dM() || j >= $self->dN())
      throw std::range_error("Index out of range");
    $self->operator[](i)[j] = val;
  }
}

%define MATSYM_EXTEND(M)
%template(MatSym ## M) nla3d::math::MatSym<M>;
%extend nla3d::math::MatSym<M> {
  %pythoncode { 
    def __str__(self):
      _list = []
      for i in range(self.dM()):
        __list = []
        for j in range(self.dN()):
          __list.append(self[i, j])
        _list.append(__list)
      return str(_list)

    def __setitem__(self, pos, val):
      i, j = pos
      self._setitem(i, j, val)

    def __getitem__(self, pos):
      i, j = pos
      return self._getitem(i, j)
  } 

  double _getitem(uint16 i, uint16 j) {
    if (i >= $self->dM() || j >= $self->dN())
      throw std::range_error("Index out of range");
    return $self->comp(i, j);
  }

  void _setitem(uint16 i, uint16 j, double val) {
    if (i >= $self->dM() || j >= $self->dN())
      throw std::range_error("Index out of range");
    $self->comp(i, j) = val;
  }
}
%enddef

MATSYM_EXTEND(3)

%rename(_addNode) nla3d::FEStorage::addNode;
%rename(_addElement) nla3d::FEStorage::addElement;
%include "FEStorage.h"
%extend nla3d::FEStorage {
  %pythoncode {
    def addNode(self, node):
        node.thisown = 0 # pass ownership to nla3d::FEStorage
        return self._addNode(node)
    def addElement(self, el):
        el.thisown = 0 # pass ownership to nla3d::FEStorage
        return self._addElement(el)
  }

  nla3d::ElementTETRA0& getElementTETRA0(uint32 _en) {
    return $self->getElement<nla3d::ElementTETRA0>(_en);
  }

  nla3d::ElementTRUSS3& getElementTRUSS3(uint32 _en) {
    return $self->getElement<nla3d::ElementTRUSS3>(_en);
  }

  nla3d::ElementQUADTH& getElementQUADTH(uint32 _en) {
    return $self->getElement<nla3d::ElementQUADTH>(_en);
  }

  nla3d::SurfaceLINETH& getElementSurfaceLINETH(uint32 _en) {
    return $self->getElement<nla3d::SurfaceLINETH>(_en);
  }
}
%include "PostProcessor.h"

%ignore nla3d::VtkProcessor::registerResults;
%include "VtkProcessor.h"

%rename(_addPostProcessor) nla3d::FESolver::addPostProcessor;
%include "FESolver.h"
%extend nla3d::FESolver {
  %pythoncode {
    def addPostProcessor(self, postprocessor):
      postprocessor.thisown = 0 # pass ownership to nla3d::FESolver
      return self._addPostProcessor(postprocessor)
  }
}

%ignore "UNDEF";
%ignore "LAST";
%include "elements/query.h"

%ignore nla3d::Element::getNodeNumber;
%ignore nla3d::Element::assembleK;
%ignore nla3d::Element::assembleM;
%ignore nla3d::Element::assembleC;
%ignore nla3d::Element::pre;
%ignore nla3d::Element::build;
%ignore nla3d::Element::post;
%ignore nla3d::ElementShape;
%include "elements/element.h"
%rename("%s") nla3d::Element::getNodeNumber;
%extend nla3d::Element {
    uint32 getNodeNumber(uint16 num) {
      if (num >= $self->getNNodes())
        throw std::range_error("Index out of range");
      return $self->getNodeNumber(num);
    }
    void setNodeNumber(uint16 num, uint32 node) {
      if (num >= $self->getNNodes())
        throw std::range_error("Index out of range");
      $self->getNodeNumber(num) = node;
    }
}

%include "Node.h"
%include "Dof.h"
%include "Mpc.h"

// use only MeshData class and readCdbFile func
%rename($ignore, %$isfunction) "";
%rename($ignore, %$isclass) "";
%rename("%s") MeshData;
%rename("%s") MeshData::getDegeneratedCells;
%rename("%s") MeshData::getCellsByAttribute;
%rename("%s") MeshData::clear;
%rename("%s") MeshData::compressNumbers;
%rename("%s") readCdbFile;
%include "FEReaders.h"
%rename("%s", %$isfunction) "";
%rename("%s", %$isclass) "";
namespace std {
   %template(vectori) vector<uint32>;
   %template(vectord) vector<double>;
   %template(vectors) vector<std::string>;
   %template(vectorv3) vector<nla3d::math::Vec<3>>;
   %template(vectorvectori) std::vector<std::vector<uint32>>;
   %template(vectorloadbc) std::vector<nla3d::loadBC>;
   %template(vectorfixbc) std::vector<nla3d::fixBC>;
   %template(vectormpcp) std::vector<Mpc*>;

   %template(map_s_fecomp) std::map<std::string, FEComponent>;
   %template(map_s_vectord) std::map<std::string, std::vector<double>>;
   %template(map_s_vectori) std::map<std::string, std::vector<uint32>>;
   %template(map_s_vectors) std::map<std::string, std::vector<std::string>>;
};

%include "elements/TRUSS3.h"
%include "elements/TETRA0.h"
%include "elements/isoparametric.h"
%include "elements/QUADTH.h"


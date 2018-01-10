// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "math/Vec.h"

namespace nla3d {
namespace math {


dVec::dVec() {

}


dVec::dVec(dVec const& rhs) {
  if(rhs.isInit()) {
    double* ptr = new double[rhs.size()];
	memcpy(ptr, rhs.data, sizeof(double) * rhs.size());
    data = ptr;
    memory_owner = true;
    _size = rhs.size();
  }
}


dVec::dVec(dVec&& rhs)
: data(rhs.data)
, memory_owner(rhs.memory_owner)
, _size(rhs.size())
{
  rhs._size = 0;
  rhs.data = nullptr;
  rhs.memory_owner = false;
}


dVec::dVec(uint32 _n, double _val) {
  reinit(_n, _val);
}


dVec::dVec(dVec& _ref, uint32 _start, uint32 _size) {
  reinit(_ref, _start, _size);
}


void dVec::reinit(uint32 _n, double _val) {
  clear();
  data = new double[_n];
  memory_owner = true;
  _size = _n;
  fill(_val);
}

void dVec::reinit(dVec& _ref, uint32 _start, uint32 _size) {
  clear();

  assert(_ref.data);
  assert(_start + _size <= _ref.size());
  data = _ref.data + _start;
  this->_size = _size;
  memory_owner = false;
}


dVec::~dVec() {
  clear();  
}

uint32 dVec::size() const {
  return _size;
}


void dVec::zero() {
  fill(0.0);
}


double* dVec::ptr() {
  assert(data);
  return data;
}


void dVec::clear() {
  if (memory_owner && data) {
    delete[] data;
  }
  memory_owner = false;
  data = nullptr;
  _size = 0;
}


void dVec::fill(double val) {
  assert(data);
  std::fill_n(data, _size, val);

}


bool dVec::isInit() const {
  return (data != nullptr);
}


double& dVec::operator[](uint32 _n) {
  assert(data);
  assert(_n < _size);
  return data[_n];
}


double dVec::operator[](uint32 _n) const {
  assert(data);
  assert(_n < _size);
  return data[_n];
}


dVec dVec::operator-() {
  assert(data);
  dVec p(size());
	for (uint32 i = 0; i < size(); i++) {
    p[i] = -data[i];
  }
  return p;
}


dVec dVec::operator+(const dVec& op) {
  assert(data);
  assert(op.data);
  assert(size() == op.size());

  dVec p(size());

	for (uint32 i = 0; i < size(); i++) {
    p[i] = data[i] + op.data[i];
  }
  return p;

}


dVec dVec::operator-(const dVec& op) {
  assert(data);
  assert(op.data);
  assert(size() == op.size());

  dVec p(size());

	for (uint32 i = 0; i < size(); i++) {
    p[i] = data[i] - op.data[i];
  }
  return p;
}


dVec dVec::operator*(const double op) {
  assert(data);

  dVec p(size());

	for (uint32 i = 0; i < size(); i++) {
    p[i] = data[i] * op;
  }
  return p;
}

dVec operator* (const double op1, const dVec& op2) {
  assert(op2.data);

  dVec p(op2.size());

	for (uint32 i = 0; i < op2.size(); i++) {
    p[i] = op2.data[i] * op1;
  }
  return p;
}


dVec dVec::operator/(const double op) {
  assert(data);

  dVec p(size());

	for (uint32 i = 0; i < size(); i++) {
    p[i] = data[i] / op;
  }
  return p;

}


dVec& dVec::operator+=(const dVec& op) {
  assert(data);
  assert(op.data);
  assert(size() == op.size());

	for (uint32 i = 0; i < size(); i++) {
    data[i] += op.data[i];
  }
  return *this;
}


dVec& dVec::operator-=(const dVec& op) {
  assert(data);
  assert(op.data);
  assert(size() == op.size());

	for (uint32 i = 0; i < size(); i++) {
    data[i] -= op.data[i];
  }
  return *this;
}

dVec& dVec::operator=(const dVec& op) {
  if (this == &op) {
    return *this;
  } else if (!op.isInit()) {
    clear();
  } else {
      if (_size != op.size())
        reinit(op.size());
      memcpy(data, op.data, sizeof(double) * op.size());
  }
  return *this;
}


bool dVec::compare(const dVec& op2, double th) {
  dVec& op1 = *this;
  assert(op1.data && op2.data);

  // first round. compare lengths
  if (op1.size() != op2.size()) {
    return false;
  }

  // second round. compare values
  for (uint32 i = 0; i < op1.size(); i++) {
    if (fabs(op1.data[i] - op2.data[i]) > th) {
      return false;
    }
  }

  return true;
}


void dVec::writeTextFormat(std::ostream& out) {
  out << size() << std::endl;
  for (uint32 i = 0; i < size(); i++) {
    out << data[i] << std::endl;
  }
}


void dVec::readTextFormat(std::istream& in) {
  uint32 _len;
  in >> _len;
  reinit(_len);
  for (uint32 i = 0; i < _len; i++) {
    in >> data[i];
  }
}

} // namespace math
} // namespace nla3d

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include <iostream>
#include <string>

namespace nla3d {
namespace math {


template<uint16 dim> 
class Vec
{
public:
	Vec() {
		assert(dim);
		memset(data,0,sizeof(double)*dim); //не очень красиво
	}
	Vec(double first, ...) {
		assert(dim);
		va_list argp;
		va_start(argp, first);
		data[0]=first;
		for (uint16 i=1; i < dim; i++) data[i]=va_arg(argp, double);
		va_end(argp);
	}
	~Vec(void) { }
  // TODO: this should be size_t instead of uint16 
	double& operator[] (uint16 n);
	const double operator[] (uint16 n) const;
	void display ();
	void zero() {
		memset(data,0,sizeof(double)*dim);
	}
	double length ();
	double qlength ();
	Vec operator+ (const Vec<dim> &op);
	Vec& operator+= (const Vec<dim> &op);
	Vec operator- ();
	Vec& operator= (const Vec<dim> &op);
	Vec operator- (const Vec<dim> &op);
	Vec operator* (const double op);
    std::string toString ();
	bool compare (Vec<dim>& V, double eps = 0.00005);
	void simple_read(std::istream& st);
	double operator* (const Vec<dim> &op);
    uint16 size() const {
      return dim;
    }
  template <uint16 dim1>
    friend std::ostream& operator<< (std::ostream& stream,const Vec<dim1>& obj);
  template <uint16 dim1>
    friend Vec operator* (const double op1, const Vec<dim1> &op2);
	double* ptr ();
private:
	double data[dim];
};
//------operator[]------------------------------------------------------
template<uint16 dim>
double& Vec<dim>::operator[] (uint16 n) {
	assert(n < dim);
	return data[n];
}
//-------operator[]-const-----------------------------------------------
template<uint16 dim>
const double Vec<dim>::operator[] (uint16 n) const {
	assert(n < dim);
	return data[n];
}
//--------display----------------------------------------------------
template<uint16 dim>
void Vec<dim>::display () {
	for (uint16 i = 0; i < dim; i++) {
		std::cout << data[i];
		if (i < dim-1) std::cout << ", ";
	}
}
//-------operator+-----------------------------------------------------
template<uint16 dim>
Vec<dim> Vec<dim>::operator+(const Vec<dim> &op) {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=data[i]+op.data[i];
	return p;
}
//-------operator+=-----------------------------------------------------
template<uint16 dim>
Vec<dim>& Vec<dim>::operator+= (const Vec<dim> &op) {
	for (uint16 i=0; i < dim; i++) this->data[i]+=op.data[i];
	return *this;
}
//-------operator- ----------------------------------------------------
template<uint16 dim>
Vec<dim> Vec<dim>::operator-() {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=-data[i];
	return p;
}
//------operator- -----------------------------------------------------
template<uint16 dim>
Vec<dim> Vec<dim>::operator-(const Vec<dim> &op) {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=data[i]-op.data[i];
	return p;
}
//------operator*------------------------------------------------------
template<uint16 dim>
Vec<dim> Vec<dim>::operator*(const double op) {
	Vec<dim> p;
	for (uint16 i=0; i < dim; i++) p[i]=data[i]*op;
	return p;
}
//------operator*------------------------------------------------------
template<uint16 dim>
double Vec<dim>::operator*(const Vec<dim> &op) {
	double p=0;
	for (uint16 i=0; i < dim; i++) p+=data[i]*op.data[i];
	return p;
}
//---------qlength----------------------------------------------------
template<uint16 dim>
double Vec<dim>::qlength() {
	double p=0;
	for (uint16 i=0; i < dim; i++) p+=data[i]*data[i];
	return p;
}
//---------length----------------------------------------------------
template<uint16 dim>
double Vec<dim>::length() {
	return sqrt(qlength());
}
//---------operator*----------------------------------------------------
template<uint16 dim1> Vec<dim1> operator*(const double op1, const Vec<dim1> &op2) {
	Vec<dim1> p;
	for (uint16 i=0; i < dim1; i++) p[i]=op2.data[i]*op1;
	return p;
}
//--------operator=------------------------------------------------------
template<uint16 dim> Vec<dim>& Vec<dim>::operator= (const Vec<dim> &op)
{
	memcpy(this->data, op.data, sizeof(double)*dim);
	return *this;
}
//---------operator<<----------------------------------------------------
template<uint16 dim1> std::ostream &operator<<(std::ostream &stream, const Vec<dim1> &obj) {
    stream << obj.data[0];
	for (uint16 i = 1; i < dim1; i++) {
		stream << " " << obj.data[i];
	}
	return stream;
}
//-------------------------------------------------------------
template<uint16 dim> std::string Vec<dim>::toString ()
{
    std::string p;
	char buff[100];
	for (uint16 i = 0; i < dim; i++) {
		sprintf_s(buff,100,"%8.5e",this->data[i]);
		p+=buff;
		if (i < dim-1) p+= ", ";
	}
	return p;
}
//-------------------------------------------------------
template<uint16 dim> double* Vec<dim>::ptr ()
{
	return data;
}
//-------------------------------------------------------
template<uint16 dim> bool Vec<dim>::compare (Vec<dim>& V, double eps)
{
	double *Dp  = (double*) data;
	double *Vp = V.ptr();
	for (uint16 j=0;j<dim;j++)
	{
		if (fabs((double)*Dp-*Vp) > eps) {
			return false;
		}
		Dp++;
		Vp++;
	}
	return true;
}
//-------------------------------------------------------
template<uint16 dim> void Vec<dim>::simple_read(std::istream& st)
{
	double *Dp = data;
	for (uint16 j=0;j<dim;j++)
	{
		st >> *Dp;
		Dp++;
	}
}



// dVec class represents a dynamically allocated double arrays with math operations support
// dVec has to modes: 1. dVec owns allocated memory; 2. dVec just point to double array memory
// allocated by another dVec. memory_owner variable is used to distinct these modes
class dVec {
  public:
    dVec();
    dVec(dVec const& rhs);
    dVec(dVec&& rhs);
    dVec(uint32 _n, double _val = 0.0);
    dVec(dVec& _ref, uint32 _start, uint32 _size);
    void reinit(uint32 _n, double _val = 0.0);
    void reinit(dVec& _ref, uint32 _start, uint32 _size);
    ~dVec();

    uint32 size() const;
    void zero();
    double* ptr();
    void clear();
    void fill(double val);

    bool isInit() const;

    double& operator[](uint32 _n);
    double operator[](uint32 _n) const;

    dVec operator-();

    dVec operator+(const dVec& op);
    dVec operator-(const dVec& op);
    dVec operator*(const double op);
    friend dVec operator* (const double op1, const dVec& op2);
    dVec operator/(const double op);

    dVec& operator+=(const dVec& op);
    dVec& operator-=(const dVec& op);

    dVec& operator=(const dVec& op);

    // for debug purpose:
    bool compare(const dVec& op2, double th = 1.0e-10);
    // write and read to/from simple text format:
    // n
    // value1
    // value2
    // ...
    // valuen
    void writeTextFormat(std::ostream& out);
    void readTextFormat(std::istream& in);

  private:

    bool memory_owner = false;
    double* data = nullptr;
    uint32 _size = 0;
};


} // namespace math
} // namespace nla3d

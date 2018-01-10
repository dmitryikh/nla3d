// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "sys.h"
#include <iostream>
#include "math/Vec.h"

#ifdef NLA3D_USE_BLAS
  #include <mkl.h>
#endif

namespace nla3d {
namespace math {

template<uint16 dimM, uint16 dimN>
class Mat
{
public:
	Mat() {
		assert(dimM && dimN);
	}
	Mat(double first, ...) {
		assert(dimM && dimN);
		va_list argp;
		va_start(argp, first);
		for (uint16 i=0; i < dimM; i++)
			for (uint16 j=0; j < dimN; j++)
			{
				if (i == 0 && j == 0) data[i][j]=first;
				else data[i][j]=va_arg(argp, double);
			}
		va_end(argp);
	}
	~Mat() { }

	Vec<dimN>& operator[] (uint16 n);
	const Vec<dimN>& operator[] (uint16 n) const;
	Mat<dimN,dimM> transpose ();
	void display ();
	Mat operator+ (const Mat &op);
	Mat& operator+= (const Mat &op);
	Mat operator- () const;
	Mat& operator= (const Mat &op);
	Mat<dimM,1>& operator= (const Vec<dimM> &op); //TODO: CHECK this operation
	Mat operator- (const Mat &op) const;
	Mat operator* (const double op);
	void Identity ();
	void zero ();
	double det();
	Vec<dimM> eigenvalues();
	Mat<dimM,dimN> inv(double det);
	Mat<dimM-1,dimN-1> cross_cut (uint16 cuti, uint16 cutj);
    std::string toString ();
	double* ptr ();
	bool compare (const Mat<dimM,dimN> &B, double eps = 1.0e-5);
	void simple_read (std::istream &st);
	//friend функции
	template <uint16 dimM1, uint16 dimN1> 
	friend std::ostream& operator<< (std::ostream& stream, const Mat<dimM1,dimN1> &obj);
	template <uint16 dimM1, uint16 dimN1, uint16 dimM2, uint16 dimN2> 
    friend Mat<dimM1,dimN2> operator* (const Mat<dimM1,dimN1> &op1, const Mat<dimM2,dimN2> &op2);
	template <uint16 dimM1, uint16 dimN1, uint16 dimM2> 
    friend Vec<dimM1> operator* (const Mat<dimM1,dimN1> &op1, const Vec<dimM2> &op2);
    // template<uint16 dimM1, uint16 dimN1>
    // friend bool matCompare (const Mat<dimM1, dimN>& mat1, const Mat<dimM1, dimN1>& mat2, const double eps);

    uint16 dM() const {
      return dimM;
    }

    uint16 dN() const {
      return dimN;
    }

	Vec<dimN> data[dimM];
private:
    //data was moved from here to public
};

//---------operator<<----------------------------------------------------------
template<uint16 dimM, uint16 dimN> std::ostream &operator<<(std::ostream &stream, const Mat<dimM,dimN> &obj) {
	for (uint16 i = 0; i < dimM; i++) {
		stream << "[" << obj.data[i] << "]" << std::endl;
	}
	return stream;
}

//----------operator[]---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
inline Vec<dimN>& Mat<dimM,dimN>::operator[] (uint16 n) {
		assert(n < dimM);
		return data[n];
}
template<uint16 dimM, uint16 dimN>
inline const Vec<dimN>& Mat<dimM,dimN>::operator[] (uint16 n) const {
	assert(n < dimM);
	return data[n];
}
//----------Identity ()---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
void Mat<dimM,dimN>::Identity () {
	assert (dimM==dimN);
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			if (i==j)
				data[i][j]=1.0f;
			else
				data[i][j]=0.0f;
}
//----------zero()---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
void Mat<dimM,dimN>::zero() {
	memset(data,0,sizeof(double)*dimM*dimN);
}
//-------------display()------------------------------------------------------
template<uint16 dimM, uint16 dimN>
void Mat<dimM,dimN>::display() {
	for (unsigned int i = 0; i < dimM; i++) {
		std::cout << "[";
		data[i].display();
		std::cout << "]" << std::endl;
	} 
}
//----------transpose()-------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimN,dimM> Mat<dimM,dimN>::transpose () {
	Mat<dimN,dimM> p;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			p[j][i]=this->data[i][j];
	return p;
}
//----------operator+(Mat)---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN> Mat<dimM,dimN>::operator+ (const Mat<dimM,dimN> &op) {
	Mat<dimM,dimN> p;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			p[i][j]=this->data[i][j]+op.data[i][j];
	return p;
}
//----------operator-()---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN> Mat<dimM,dimN>::operator- () const {
	Mat<dimM,dimN> p;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			p[i][j] = -this->data[i][j];
	return p;
}
//----------operator-(Mat)---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN> Mat<dimM,dimN>::operator- (const Mat<dimM,dimN> &op) const {
	Mat<dimM,dimN> p;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			p[i][j]=this->data[i][j]-op.data[i][j];
	return p;
}
//-----------operator*(double)----------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN> Mat<dimM,dimN>::operator*(const double op)
{
	Mat<dimM,dimN> p;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			p[i][j]=this->data[i][j]*op;
	return p;
}
//----------operator+=---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN>& Mat<dimM,dimN>::operator+= (const Mat<dimM,dimN> &op) {
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			this->data[i][j]+=op.data[i][j];
	return *this;
}
//----------operator=(Mat)---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN>& Mat<dimM,dimN>::operator= (const Mat<dimM,dimN> &op)
{
	for (uint16 i = 0; i < dimM; i++)
		this->data[i]=op.data[i];
	return *this;
}
//----------operator=(Vec)---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,1>& Mat<dimM,dimN>::operator= (const Vec<dimM> &op)
{
	assert(dimN == 1); //только для матрицы-столбца
	for (uint16 i = 0; i < dimM; i++)
		this->data[i][0] = op[i];
	return *this;
}
//----------operator*(Mat)---------------------------------------------------------
template<uint16 dimM1, uint16 dimN1, uint16 dimM2, uint16 dimN2> Mat<dimM1,dimN2> operator*(const Mat<dimM1,dimN1> &op1, const Mat<dimM2,dimN2> &op2) {
	assert(dimN1 == dimM2);
	Mat<dimM1, dimN2> p;
	double element;
	double *ptr1 = (double*) op1.data;
	double *ptr2 = (double*) op2.data;
	for (uint16 i=0; i < dimM1; i++)
		for (uint16 j=0; j < dimN2; j++)
		{
			element=0.0f;
			for (uint16 l=0; l < dimN1; l++) element+=ptr1[i*dimN1+l]*ptr2[l*dimN2+j];
				//(op1.data[i])[l]*(op2.data[l])[j];
			p[i][j]=element;
		}
	return p;


	//assert(dimN == dimM2);
	//Mat<dimM, dimN2> p;
	//double element;
	//for (uint16 i=0; i < dimM; i++)
	//	for (uint16 j=0; j < dimN2; j++)
	//	{
	//		element=0.0f;
	//		for (uint16 l=0; l < dimN; l++) element+=(op1.data[i])[l]*(op2.data[l])[j];
	//		p[i][j]=element;
	//	}
	//return p;
}
//-----------operator*(Vec)---------------------------------------------------------
template <uint16 dimM1, uint16 dimN1, uint16 dimM2> Vec<dimM1> operator* (const Mat<dimM1,dimN1> &op1, const Vec<dimM2> &op2) {
	assert(dimN1==dimM2);
	Vec<dimM1> p;
	double el;
	for (uint16 i=0; i < dimM1; i++)
	{
		el = 0.0f;
		for (uint16 j=0; j < dimN1; j++)
			el += op1.data[i][j]*op2[j];
		p[i] = el;
	}
	return p;
}
//--------------------------------------------------------------------------------
template<uint16 dimM, uint16 dimN>
std::string Mat<dimM,dimN>::toString()
{
    std::string p;
	for (unsigned int i = 0; i < dimM; i++) {
		p+= "[";
		p+=data[i].toString();
		p+= "]";
		if (i!=dimM-1) p+= "\n";
	} 
	return p;
}
//-------------------------------------------------------------------------------
template<> double Mat<1,1>::det();
template<> double Mat<2,2>::det();
template<> double Mat<3,3>::det();

template<uint16 dimM, uint16 dimN>
double Mat<dimM,dimN>::det()
{
	//разложение по строке i=0
	assert(dimM == dimN);
	double det = 0.0;
	double koef = 1;
	for (uint16 j=0; j < dimN; j++)
	{
		Mat<dimM-1,dimM-1> tmp;
		det += koef*data[0][j]*cross_cut(0,j).det();
		koef *= -1;
	}
	return det;
}

template<> Mat<1, 1> Mat<1,1>::inv(double det);
template<> Mat<2, 2> Mat<2,2>::inv(double det);
template<> Mat<3, 3> Mat<3,3>::inv(double det);

template<uint16 dimM, uint16 dimN>
Mat<dimM, dimN> Mat<dimM,dimN>::inv(double det)
{
	//через матрицу алг. дополнений
	assert(dimM == dimN);
	Mat<dimM,dimN> C;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			C[i][j] = npow(-1,i+j)*cross_cut(i,j).det();  //CHECK
	C = C.transpose()*(1/det);
	return C;
}

template<uint16 dimM, uint16 dimN>
Mat<dimM-1,dimN-1> Mat<dimM,dimN>::cross_cut (uint16 cuti, uint16 cutj)
{
	Mat<dimM-1, dimN-1> tmp;
	for (uint16 ii=0; ii < dimM-1; ii++)
		for (uint16 jj=0; jj < dimN-1; jj++)
			tmp[ii][jj] = data[(ii>=cuti)?ii+1:ii][(jj>=cutj)?jj+1:jj];
	return tmp;
}
template<uint16 dimM, uint16 dimN>
Vec<dimM>  Mat<dimM,dimN>::eigenvalues()
{
	assert(dimM == 3 && dimN == 3); //TODO: пока только для матриц 3 на 3
	double I1 = data[0][0] + data[1][1] + data[2][2];
	double I2 = data[0][0]*data[1][1] + data[1][1]*data[2][2] + data[0][0]*data[2][2] - 
		data[0][1]*data[0][1] - data[0][2]*data[0][2] - data[1][2]*data[1][2];
	double I3 = data[0][0]*data[1][1]*data[2][2] + 2*data[0][1]*data[1][2]*data[0][2] - 
		data[0][0]*data[1][2]*data[1][2]-data[1][1]*data[0][2]*data[0][2]-data[2][2]*data[0][1]*data[0][1];

	double a = -I1;
	double b = I2;
	double c = -I3;
	 // код из http://ru.wikipedia.org/wiki/%D0%A2%D1%80%D0%B8%D0%B3%D0%BE%D0%BD%D0%BE%D0%BC%D0%B5%D1%82%D1%80%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B0%D1%8F_%D1%84%D0%BE%D1%80%D0%BC%D1%83%D0%BB%D0%B0_%D0%92%D0%B8%D0%B5%D1%82%D0%B0
	// для трех вещественных корней
	// x*x*x + a * x * x + b * x  + c == 0
	double p = b - a * a / 3;
	double q = 2 * a * a * a / 27 - a * b / 3 + c;
	double A = sqrt(- 4 * p / 3); 
 
	double c3phi = - 4 * q / (A * A * A);
 
	double phi = acos(c3phi) / 3;  
 
	double root1 = A * cos(phi) - a / 3;
	double root2 = A * cos(phi + 2 * M_PI / 3) - a / 3;
	double root3 = A * cos(phi - 2 * M_PI / 3) - a / 3;

	//сортируем
	double s1,s2,s3;
	if (root1 > root2 && root1 > root3)
	{
		s1 = root1;
		if (root2 > root3)
		{
			s2 = root2;
			s3 = root3;
		}
		else
		{
			s2 = root3;
			s3 = root2;
		}
	}
	else if (root2 > root1 && root2 > root3)
	{
		s1 = root2;
		if (root1 > root3)
		{
			s2 = root1;
			s3 = root3;
		}
		else
		{
			s2 = root3;
			s3 = root1;
		}
	}
	else
	{
		s1 = root3;
		if (root1 > root2)
		{
			s2 = root1;
			s3 = root2;
		}
		else
		{
			s2 = root2;
			s3 = root1;
		}
	}
	return Vec<3>(s1,s2,s3);
}

template<uint16 dimM, uint16 dimN>
double* Mat<dimM,dimN>::ptr ()
{
	return data[0].ptr();
}


template<uint16 dimM, uint16 dimN>
bool Mat<dimM,dimN>::compare (const Mat<dimM,dimN> &B, double eps) {
	double *Dp = (double*) data;
	// double *Bp = B.ptr();
    double *Bp = (double*) B.data;
	for (uint16 i = 0; i < dimM; i++) {
		for (uint16 j = 0; j < dimN; j++) {
			if (fabs(*Dp-*Bp) > eps) {
        LOG(INFO) << "Mat[" << i << "][" << j << "]: " << *Dp << " != " << *Bp;
				return false;
			}
			Dp++;
			Bp++;
		}
	}
	return true;
}


template<uint16 dimM, uint16 dimN>
void Mat<dimM,dimN>::simple_read (std::istream &st)
{
	double *Bp = (double*) data;
	for (uint16 i=0;i<dimM;i++)
	{
		for (uint16 j=0;j<dimN;j++)
		{
			st >> *Bp;
			Bp++;
		}
	}
}


class dMat;
class dMat_interface
{
public:
	dMat_interface ():ptr(NULL),M(0),N(0),row(0) 
	{	}
	double& operator[] (uint16 col)
	{
		assert(ptr);
		assert (col < N);
		assert (row < M);
		return ptr[row*N + col];
	}
	friend class dMat;
private:
	uint16 M, N;
	uint32 row;
	double *ptr;
};

//dynamic matrix to pass arbitrary matrix to functions as arguments
class dMat
{
public:
	dMat(uint16 dim_m, uint16 dim_n) : dimM(0),dimN(0),data(NULL)
	{
		if (dim_m && dim_n)
			resize(dim_m, dim_n);
	}
	dMat(uint16 dim_m, uint16 dim_n, double first, ...) : dimM(0),dimN(0),data(NULL) {
    va_list argp;
		if (dim_m && dim_n)
		{
			resize(dim_m, dim_n);
			va_start(argp, first);
			data[0] = first;
			for (uint16 i=1; i < dimM*dimN; i++)
				data[i]=va_arg(argp, double);
			va_end(argp);
		}
	}
	dMat(const dMat &from) : dimM(0),dimN(0),data(NULL)
	{
		operator=(from);
	}
	dMat_interface& operator[] (uint16 row) {
        //TODO: it seems unsafe for parallel read access to dMat..
		dmat_int.row = row;
		return dmat_int;
	}
	void resize (uint16 dim_m, uint16 dim_n)
	{
		if (data) delete[] data;
		data = new double[dim_m*dim_n];
		dmat_int.ptr = data;
		dmat_int.M = dim_m;
		dmat_int.N = dim_n;
		dimM = dim_m;
		dimN = dim_n;
		zero();
	}

	void fill (double first, ...) {
		va_list argp;
		va_start(argp, first);
		data[0] = first;
		for (uint16 i=1; i < dimM*dimN; i++)
			data[i]=va_arg(argp, double);
		va_end(argp);
	}
	~dMat()
	{
		if (data) delete[] data;
		data = NULL;
	}

	void zero ()
	{
		memset(data,0,sizeof(double)*dimM*dimN);
	}
	template <uint16 M, uint16 N> Mat<M,N> toMat();
	template <uint16 M> Vec<M> toVec(uint16 col=0);
	template <uint16 M, uint16 N> void cpMat(Mat<M,N> &mat);
	template <uint16 M> void cpVec(Vec<M> &vec, uint16 col=0);
	dMat& operator= (const dMat &from)
	{
		resize(from.dimM, from.dimN);
		memcpy(this->data, from.data, sizeof(double)*dimM*dimN);
		return *this;
	}

	uint16 dM ()
	{
		return dimM;
	}
	uint16 dN ()
	{
		return dimN;
	}

    double* ptr()
    {
        return data;
    }

	friend std::ostream& operator<< (std::ostream& stream, dMat &obj);
private:
	uint16 dimM;
	uint16 dimN;
	double *data;
	dMat_interface dmat_int;
};


//копирует dMat в правый верхний угол Mat<M,N>
template <uint16 M, uint16 N>
Mat<M,N> dMat::toMat()
{
	assert(M >= dimM && N >= dimN); //Mat >= dMat
	Mat<M,N> tmp;
	for (uint16 i = 0; i < M; i++)
		for (uint16 j = 0; j < N; j++)
			tmp[i][j] = (*this)[i][j];
	return tmp;
}

template <uint16 M> 
Vec<M> dMat::toVec(uint16 col)
{
	assert(M >= dimM && dimN > col); //Vec >= dMat по M 
	Vec<M> tmp;
	for (uint16 i = 0; i < M; i++)
		tmp[i] = (*this)[i][col];
	return tmp;
}

template <uint16 M, uint16 N> 
void dMat::cpMat(Mat<M,N> &mat)
{
	assert(dimM >= M && dimN >= N); //dMat >= Mat
	for (uint16 i = 0; i < M; i++)
		for (uint16 j = 0; j < N; j++)
			(*this)[i][j] = mat[i][j];
}

//функция копирует вектор в столбец col матрицы dMat
// col - от нуля
template <uint16 M> 
void dMat::cpVec(Vec<M> &vec, uint16 col)
{
	assert(dimM >= M && dimN > col);
	for (uint16 i = 0; i < M; i++)
		(*this)[i][col] = vec[i];
}


//-------------------------------------------------------------
//	MatSym
//-------------------------------------------------------------
template<uint16 dimM>
class MatSym {
public:
	MatSym () {
		assert(dimM);
	}
	double* ptr() {
		return data;
	}

	const double* ptr() const {
		return data;
	}

	void zero() {
		memset((void*) data, 0, sizeof(double)*getLength());
	}
	
  // indexing from 0
	double& comp (uint16 i, uint16 j) {
		uint16 b;
		if (i > j) {
			b = j;
			j = i;
			i = b;
		}
		b = (2*dimM-(i-1))/2*i + (j-i);
		b = dimM*i - (i-1)*i/2 + (j-i);
		return data[b];
	}
	
	uint16 getLength () {
		return dimM*(dimM+1)/2;
	}
  Mat<dimM, dimM> toMat();
	void simple_read (std::istream &st);
	bool compare (MatSym<dimM> &B, double eps = 1.0e-5);
	MatSym& operator+= (const MatSym &op);
    uint16 dM() const {
      return dimM;
    }
    uint16 dN() const {
      return dimM;
    }

	double data[dimM*(dimM+1)/2];
};

template<uint16 dimM>
void MatSym<dimM>::simple_read (std::istream &st)
{
	double *Bp = (double*) data;
	uint16 l = getLength();
	for (uint16 i=0;i<l;i++)
	{
			st >> *Bp;
			Bp++;
	}
}

template<uint16 dimM>
bool MatSym<dimM>::compare (MatSym<dimM> &B, double eps)
{
	double *Dp = (double*) data;
	double *Bp = B.ptr();
	uint16 l = getLength();
	for (uint16 i=0;i<l;i++)
	{
		if (fabs(*Dp-*Bp) > eps) {
			return false;
		}
		Dp++;
		Bp++;
	}
	return true;
}

template<uint16 dimM>
Mat<dimM, dimM> MatSym<dimM>::toMat() {
  uint16 ind = 0;
  Mat<dimM,dimM> mat;
  for (uint16 i = 0; i < dimM; i++) {
    for (uint16 j = i; j < dimM; j++) {
      mat[i][j] = data[ind];
      mat[j][i] = data[ind];
      ind++;
    }
  }
  return mat;
}


template<uint16 dimM>
MatSym<dimM>& MatSym<dimM>::operator+= (const MatSym<dimM> &op) {
  double* thisp = this->ptr();
  const double* opp = op.ptr();
  for (uint16 i = 0; i < getLength(); i++) {
    *thisp += *opp;
    thisp++;
    opp++;
  }

	return *this;
}

// [R] = coef * [B]^T*[D]*[B]
// coef - double scalar
// [B] - common matrix (dimM x dimN)
// [D] - symmetric matrix (dimM x dimM)
// [R] - symmetrix matrix (dimN x dimN)
template<uint16 dimM,uint16 dimN>
void matBTDBprod (Mat<dimM,dimN> &B, MatSym<dimM> &D, double coef, MatSym<dimN> &R) 
{
	Mat<dimN,dimM> A;
	A.zero();
	double *Ap = A.ptr();
	double *Bp = B.ptr();
	double *Dp = D.ptr();
	double *Rp = R.ptr();
	uint16 i,j,k;

	//A = BT*D
		// BT*Du + BT*Dd
	for (k=0;k<dimM;k++)
	{
		j = k; // include diag elem
		while (j < dimM)
		{
			for (i=0;i<dimN;i++)
			{
				Ap[i*dimM+j] += Bp[k*dimN+i]*(*Dp);
			}
			Dp++;
			j++;
		}
	}
		// BT*DuT
		Dp = D.ptr();
	for (j=0;j<dimM;j++)
	{
		k = j + 1; Dp++; //leave the diagonal
		while (k < dimM)
		{
			for (i=0;i<dimN;i++)
			{
				Ap[i*dimM+j] += Bp[k*dimN+i]*(*Dp);
			}
			Dp++;
			k++;
		}
	}

	//R = A*B
	for (i=0;i<dimN;i++)
	{
		for (j=i;j<dimN;j++)
		{
			for (k=0;k<dimM;k++)
			{
				*Rp += Ap[i*dimM+k]*Bp[k*dimN+j]*coef;
			}
			Rp++;
		}
	}

}


template<uint16 dimM,uint16 dimN>
void matBTVprod(Mat<dimM,dimN> &B, Vec<dimM> &V, double coef, Vec<dimN> &R)
{
#ifndef NLA3D_USE_BLAS
	double *Bp = B.ptr();
	uint16 i,j;
	for (i=0;i<dimN;i++)
		for (j=0;j<dimM;j++)
			R[i] += Bp[j*dimN+i]*V[j]*coef;
#else
  cblas_dgemv(CblasRowMajor, CblasTrans, dimM, dimN, coef, B.ptr(), dimN, V.ptr(), 1, 0.0, R.ptr(), 1);
#endif
}

template<uint16 dimM,uint16 dimN>
void matBVprod(Mat<dimM,dimN> &B,Vec<dimN> &V, double coef, Vec<dimM> &R) {
#ifndef NLA3D_USE_BLAS
	double *Bp = B.ptr();
	uint16 i,j;
	for (i=0;i<dimM;i++)
		for (j=0;j<dimN;j++)
			R[i] += Bp[i*dimN+j]*V[j]*coef;
#else
  cblas_dgemv(CblasRowMajor, CblasNoTrans, dimM, dimN, coef, B.ptr(), dimN, V.ptr(), 1, 0.0, R.ptr(), 1);
#endif
}

template<uint16 dimM>
void matBVprod(MatSym<dimM> &B,Vec<dimM> &V, double coef, Vec<dimM> &R) {
  // TODO: this this inefficient version
  // convert to regular matrix
  Mat<dimM, dimM> mB = B.toMat();
	double *Bp = mB.ptr();
	uint16 i,j;
	for (i=0;i<dimM;i++)
		for (j=0;j<dimM;j++)
			R[i] += Bp[i*dimM+j]*V[j]*coef;
}


template<uint16 dimM1,uint16 dimN1,uint16 dimN2>
void matABprod(Mat<dimM1,dimN1> &A, Mat<dimN1,dimN2> &B, const double coef, Mat<dimM1,dimN2> &R) {
#ifndef NLA3D_USE_BLAS
	uint16 i,j,k;
	double *Rp = R.ptr();
	double *Ap = A.ptr();
	double *Bp = B.ptr();
	for (i=0;i<dimM1;i++) {
		for (j=0;j<dimN2;j++) {
			for (k=0;k<dimN1;k++)
				*Rp += Ap[i*dimN1+k]*Bp[k*dimN2+j]*coef;
			Rp++;
		}
	}
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimM1, dimN2, dimN1, coef, A.ptr(), dimN1, B.ptr(), dimN2, 0.0, R.ptr(), dimN2);
#endif
}

template<uint16 dimM1,uint16 dimN1,uint16 dimN2>
void matATBprod(Mat<dimM1,dimN1> &A, Mat<dimM1,dimN2> &B, const double coef, Mat<dimN1,dimN2> &R)  {
#ifndef NLA3D_USE_BLAS
	uint16 i,j,k;
	double *Rp = R.ptr();
	double *Ap = A.ptr();
	double *Bp = B.ptr();
	for (i=0;i<dimN1;i++)
	{
		for (j=0;j<dimN2;j++)
		{
			for (k=0;k<dimM1;k++)
				*Rp += Ap[k*dimN1+i]*Bp[k*dimN2+j]*coef;
			Rp++;
		}
	}
#else
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dimN1, dimN2, dimM1, coef, A.ptr(), dimN1, B.ptr(), dimN2, 0.0, R.ptr(), dimN2);
#endif
}

} // namespace math
} // namespace nla3d

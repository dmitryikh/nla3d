#pragma once
#include "sys.h"

class Mat_band_rm;
class Mat_band_cm;
// чтобы пользоваться классом Mat_Band вот так M[i][j]
// _rm - row major
class Mat_Band_interface_rm {
public:
	Mat_Band_interface_rm (uint32 nrow, uint32 nband):ptr(0),nRow(nrow),nBand(nband) 
	{  }
	double& operator[] (uint32 col)
	{
		uint32 i,j;
		assert (ptr);
		assert (col > 0);
		assert (row > 0);
		if (col >= row) //верхний треугольник
		{
			i = row-1;
			j = col-row;
		}
		else
		{
			i = col-1;
			j = row-col;
		}
		assert (i < nRow);
		assert (j < nBand);
		return ptr[i*nBand + j];
	}
	friend class Mat_Band_rm;
private:
	uint32 nBand, nRow;
	uint32 row;
	double *ptr;
};

// _cm - column major
class Mat_Band_interface_cm {
public:
	Mat_Band_interface_cm (uint32 n, uint32 nband):ptr(0),n(n),nBand(nband) 
	{  }
	double& operator[] (uint32 col)
	{
		uint32 i,j;
		assert (ptr);
		assert (col > 0);
		assert (row > 0);
		if (col >= row) //верхний треугольник
		{
			j=col-1;
			i=nBand-(col-row)-1;
		}
		else
		{
			j=row-1;
			i=nBand-(row-col)-1;
		}
		assert (i < nBand);
		assert (j < n);
		return ptr[j*nBand + i];
	}
	friend class Mat_Band_cm;
private:
	uint32 nBand, n;
	uint32 row;
	double *ptr;
};

//класс представляет механизмы и хранение для ленточной матрицы
// _rm - row major
class Mat_Band_rm {
public:
	Mat_Band_rm (uint32 n, uint32 nband);
	Mat_Band_interface_rm& operator[] (uint32 row) {
		intface.row = row;
		return intface;
	}
	~Mat_Band_rm() {
		if (ptr) delete[] ptr;
	}
	void zero () {
		memset(ptr, 0, nBand*n*sizeof(double));
	}
	double* Ptr () {
		return ptr;
	}
	friend std::ostream& operator<< (std::ostream &stream, const Mat_Band_rm &mat);
private:
	uint32 nBand, n;
	Mat_Band_interface_rm intface;
	double *ptr;
};

class Mat_Band_cm {
public:
	Mat_Band_cm (uint32 n, uint32 nband);
	Mat_Band_interface_cm& operator[] (uint32 row) {
		intface.row = row;
		return intface;
	}
	~Mat_Band_cm() {
		if (ptr) delete[] ptr;
	}
	void zero () {
		memset(ptr, 0, nBand*n*sizeof(double));
	}
	double* Ptr () {
		return ptr;
	}
	void toUTriangular (double *p);
	Mat_Band_cm& operator= (const Mat_Band_cm &op);
	friend std::ostream& operator<< (std::ostream &stream, const Mat_Band_cm &mat);
	string toString ();

private:
	uint32 nBand, n;
	Mat_Band_interface_cm intface;
	double *ptr;
};

//длинный вектор
class Vec_Long {
public:
	Vec_Long (uint32 nrow)
	{
		assert (nrow > 0);
		nRow=nrow;
		ptr = new double[nRow];
		zero();
	}
	~Vec_Long () {
		if (ptr) delete[] ptr;
	}
	double& operator[] (uint32 i) {
		assert ((i>0) && (i-1<nRow));
		return ptr[i-1];
	}
	double* Ptr () {
		return ptr;
	}
	void zero () {
		memset(ptr, 0, nRow*sizeof(double));
	}
	double fsum () {
		double s = 0.0f;
		for (uint16 i = 0; i < nRow; i++)
			s+=fabs(ptr[i]);
		return s;
	}
	Vec_Long& operator= (Vec_Long &op)
	{
		assert(this->nRow == op.nRow);
		memcpy(this->ptr, op.ptr, sizeof(double)*nRow);
		return *this;
	}
	void multiply(const double op)
	{
		for (uint16 i = 0; i < nRow; i++)
			ptr[i]*=op;
	}
	string toString() {
		string p;
		char buff[100];
		for (uint16 i = 0; i < nRow; i++) {
			sprintf_s(buff,100,"%f",this->ptr[i]);
			p+=buff;
			if (i < nRow-1) p+= ", ";
		}
		return p;
	}
	//прибавление к this вектору op вектор, такой же размерностью
	void add (Vec_Long *op) {
		assert(this->nRow == op->nRow);
		for (uint32 i = 0; i < nRow; i++) this->ptr[i]+=op->ptr[i];
	}
	friend std::ostream& operator<< (std::ostream &stream, const Vec_Long &vec);
private:
	double *ptr;
	uint32 nRow;
};


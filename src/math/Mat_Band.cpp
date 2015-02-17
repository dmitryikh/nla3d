#include "Mat_Band.h"


Mat_Band_rm::Mat_Band_rm (uint32 n, uint32 nband):intface(n, nband)
{
	assert (nband>0 && n >0);
	ptr = new double[nband*n];
	nBand=nband;
	this->n=n;
	intface.ptr=this->ptr;
	zero();
}

std::ostream& operator<< (std::ostream &stream, const Mat_Band_rm &mat) {
	for (uint32 i=0; i < mat.n; i++)
	{
		for (uint32 j=0; j < mat.nBand; j++)
			stream << mat.ptr[i*mat.nBand+j] << '\t';
		stream << std::endl;
	}
	return stream;
}

Mat_Band_cm::Mat_Band_cm (uint32 n, uint32 nband):intface(n, nband)
{
	assert (nband>0 && n >0);
	ptr = new double[nband*n];
	nBand=nband;
	this->n=n;
	intface.ptr=this->ptr;
	zero();
}

std::ostream& operator<< (std::ostream &stream, const Mat_Band_cm &mat) {
	for (uint32 i=0; i < mat.nBand; i++)
	{
		for (uint32 j=0; j < mat.n; j++)
			stream << mat.ptr[j*mat.nBand+i] << '\t';
		stream << std::endl;
	}
	return stream;
}

string Mat_Band_cm::toString ()
{
	string p;
	char buff[100];
	for (uint32 i=0; i < nBand; i++)
	{
		p+="[";
		for (uint32 j=0; j < n; j++)
		{
			sprintf_s(buff,100,"%8.5e",ptr[j*nBand+i]);
			p+= buff;
			if (j!=n-1) p+=",\t";
		}
		p+="]";
		if (i!=nBand-1) p+="\n";
	}
	return p;
}
void Mat_Band_cm::toUTriangular (double *p)
{
	//memset(p, 0, n*(n+1)/2*sizeof(double));
	for (uint32 col=0; col < n; col++)
	{
		uint32 tri_n = col+1; // сколько элементов надо скопировать в данном столбце

		uint32 n_copy = tri_n;
		if (n_copy>nBand) n_copy=nBand;

		uint32 pn = col*(1+col)/2     +     (tri_n-n_copy);
		// указатель на 1 стр.        //смещение
		uint32 ps=(col+1)*(nBand)-1-(n_copy-1);
		memcpy(p+pn, ptr+ps, n_copy*sizeof(double)); 
	}
}

Mat_Band_cm& Mat_Band_cm::operator= (const Mat_Band_cm &op)
{
	assert(this->n == op.n);
	assert(this->nBand == op.nBand);
	memcpy(this->ptr, op.ptr, n*nBand*sizeof(double));
	return *this;
}

std::ostream& operator<< (std::ostream &stream, const Vec_Long &vec) {
	for (uint32 i=0; i < vec.nRow; i++)
	{
			stream << vec.ptr[i] << '\t';
	}
	return stream;
}

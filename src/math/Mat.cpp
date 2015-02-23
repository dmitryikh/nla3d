#include "Mat.h"

//dMat dMat_tmp(1,1);

template<>
double Mat<1,1>::det()
{
	return data[0][0];
}
template<>
double Mat<2,2>::det()
{
	return data[0][0]*data[1][1]-data[0][1]*data[1][0];
}
template<>
double Mat<3,3>::det()
{
	return data[0][0]*(data[1][1]*data[2][2]-data[1][2]*data[2][1])-data[0][1]*(data[1][0]*data[2][2]-data[1][2]*data[2][0])+data[0][2]*(data[1][0]*data[2][1]-data[1][1]*data[2][0]);
}

//TODO: check det == 0
template<>
Mat<1,1> Mat<1,1>::inv(double det)
{
	return (1.0/det);
}
template<>
Mat<2,2> Mat<2,2>::inv(double det)
{
	Mat<2,2> tmp(data[1][1],-data[0][1],-data[1][0],data[0][0]);
	return tmp*(1.0/det);
}
template<>
Mat<3,3> Mat<3,3>::inv(double det)
{
	Mat<3,3> tmp(data[1][1]*data[2][2]-data[1][2]*data[2][1],data[0][2]*data[2][1]-data[0][1]*data[2][2],data[0][1]*data[1][2]-data[0][2]*data[1][1],
				data[2][0]*data[1][2]-data[1][0]*data[2][2],data[0][0]*data[2][2]-data[0][2]*data[2][0],data[1][0]*data[0][2]-data[0][0]*data[1][2],
				data[1][0]*data[2][1]-data[1][1]*data[2][0],data[0][1]*data[2][0]-data[0][0]*data[2][1],data[0][0]*data[1][1]-data[0][1]*data[1][0]);
	return tmp*(1.0/det);
}


//---------operator<<----------------------------------------------------------
std::ostream &operator<<(std::ostream &stream, dMat &obj) {
	for (uint16 i = 0; i < obj.dimM; i++)
	{
		std::cout << "[";
		for (uint16 j = 0; j < obj.dimN; j++)
			std::cout << obj[i][j] << " ";
		std::cout << "]" << endl;
	}
	return stream;
}

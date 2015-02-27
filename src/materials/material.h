#pragma once
#include <string>
#include "math\Vec.h"
#include "math\Mat.h"
#include "States.h"

class MaterialFactory;

// Material must have named material constants
class Material {
public:
	enum mat_func_deriv {
		AL_1	=	0,
		AL_2	=	1,
		AL_11	=	2,
		AL_12	=	3,
		AL_22	=	4
	};

	Material () : MC(NULL), numC(0), code(0) { 
	}
	Material (uint16 num_c);
	~Material()   // TODO: discover the virtual destructor
	{
		if (MC) delete[] MC; 
		MC = NULL;
	}

	virtual string toString();

	string getName();
	uint16 getCode ();
	double& Ci (uint16 i);
	double& getCstr (char* mname);
	uint16 getNumC ();
	void read_from_stream (istream &str);

	static double getJ(const double* C);
	static void getC_inv(const double* C, const double J, double* C_inv);
	static const double I[6];

  friend class MaterialFactory;
protected:
	void register_mat_const(uint16 num, ...);
	double* MC;
	vector<string> MC_names;
	uint16 numC;
	uint16 code;
	string name;
};


// ---=== FUNCTIONS ===--- //
inline string Material::getName()
{
	return name;
}

inline uint16 Material::getCode ()
{
	return code;
}

inline double& Material::Ci (uint16 i)
{
	assert(i < numC);
	return MC[i];
}

inline uint16 Material::getNumC ()
{
	return numC;
}


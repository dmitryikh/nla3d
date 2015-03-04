
#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <strstream>
#include <stdarg.h>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#ifdef linux
  #include <string.h>
  #include <math.h>
#endif

// MSVC doen's support C99 standard
#ifdef linux
  #define sprintf_s snprintf
#endif

using namespace std;
//#ifdef VC_COMPILER
#define int8 char //-127 to +127
#define uint8 unsigned char //0 to +255
#define int16 short //-32 767 to +32 767
#define uint16 unsigned short //0 to +65 535
#define int32 int //-2 147 483 647 to +2 147 483 647
#define uint32 unsigned int // 0 to +4 294 967 295
#define int64 long long //-9 223 372 036 854 775 807 to +9 223 372 036 854 775 807
#define uint64 unsigned long long //0 to +18 446 744 073 709 551 615

//#else
//#error Dont know integers types
//#endif

#define M_PI       3.14159265358979323846

enum el_component 
{
	COMP_UNDEF,

	E_X,
	E_Y,
	E_Z,
	E_XY,
	E_XZ,
	E_YZ,
	E_VOL,
	E_1,
	E_2,
	E_3,

	S_X,
	S_Y,
	S_Z,
	S_XY,
	S_XZ,
	S_YZ,
	S_P,
	S_1,
	S_2,
	S_3,

	U_X,
	U_Y,
	U_Z,

	POS_X,
	POS_Y,
	POS_Z,

	COMP_LAST
};

const char* const el_component_labels[]={"UNDEFINED","E_X","E_Y","E_Z","E_XY","E_XZ","E_YZ","E_VOL","E_1","E_2","E_3","S_X","S_Y","S_Z","S_XY","S_XZ","S_YZ","S_P","S_1","S_2","S_3","U_X","U_Y","U_Z","POS_X","POS_Y","POS_Z","COMP_LAST"};

enum el_tensor
{
	TENS_UNDEF,
	TENS_COUCHY,
  TENS_PK2, // second Piola-Kirchgoff stress tensor (symmetric 3x3)
  TENS_E,
  TENS_C
};

const char* const el_tensor_labels[]={"UNDEFINED","COUCHY"};

// use it in material matrix creator
#define ANALYSIS_3D 1
#define ANALYSIS_2D_PLANE_STRESS 2
#define ANALYSIS_2D_PLANE_STRAIN 3
#define ANALYSIS_2D_AXISYMMETRIC 4

#define GP_MEAN 100 //среднее значение по элементу

//ключ закрепления != 0  !
#define D_UX 1
#define D_UY 2
#define D_UZ 3

#define F_X 4
#define F_Y 5
#define F_Z 6


// TODO: need to wrap this enum into some namespace
enum tensorComponents {
	M_XX =  0,
	M_XY =  1,
	M_XZ =  2,
	M_YY =  3,
	M_YZ =  4,
	M_ZZ =  5
};
// TODO: need to wrap somewhere
const tensorComponents defaultTensorComponents[6] = {M_XX, M_XY, M_XZ, M_YY, M_YZ, M_ZZ};

const char SYS_VERSION[] = "1.1";
const char SYS_DATA[] = "30.01.12";
const char log_file_name[]="log.txt";
const bool debug_mode = true;

class FE_Storage;

class Log_opts 
{
public:
	Log_opts() 
	{
    //TODO: mutex
		//output_lock=CreateMutex(NULL, FALSE, NULL);
		ofstream file(log_file_name,ios::trunc);
		file.close();
		
	};
	//TODO: mutex 
  //HANDLE output_lock;
};

void warning(const char* logline, ...);
void debug(const char* logline, ...);
void debug(string &str);
void log(string &str);
void log(const char* logline, ...);
void error(const char* logline, ...);
void echo (const char* logline, ...);
void echolog(const char* logline, ...);
uint32 tick();

string IntToStr (uint32 dig);
int32 npow(int16 dig, uint16 power);

vector<string> read_tokens(char *input);
void del_spaces (string &str);

class Timer
{
public:
	Timer(bool _start = false) : start_time(0), end_time(0)
	{
		if (_start)
			start();
	}
	void start()
	{
		start_time = clock();
		end_time = start_time;
	}
	double stop()
	{
		end_time = clock();
		return time();
	}
	double time()
	{
		return ((double)end_time - start_time) / CLOCKS_PER_SEC;
	}
private:
	clock_t start_time;
	clock_t end_time;
};

uint16 str2dof (string dof_name);

char* getCmdOption(char ** begin, char ** end, const std::string & option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);
vector<char*> getCmdManyOptions(char ** begin, char ** end, const std::string & option); 

struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};

string getFileNameFromPath(const string filename);

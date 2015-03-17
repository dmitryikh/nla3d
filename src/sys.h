// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <strstream>
#include <stdarg.h>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#ifdef linux
  #include <string.h>
#endif

// MSVC doent's support C99 standard
#ifdef linux
  #define sprintf_s snprintf
#endif

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

namespace nla3d {

const char SYS_VERSION[] = "1.3";
const char SYS_DATA[] = "16.03.15";
const char log_file_name[] = "log.txt";
const bool debug_mode = true;

//class FEStorage;

class Log_opts 
{
public:
	Log_opts() 
	{
    //TODO: mutex
		//output_lock=CreateMutex(NULL, FALSE, NULL);
    std::ofstream file(log_file_name, std::ios::trunc);
		file.close();
		
	};
	//TODO: mutex 
  //HANDLE output_lock;
};

void warning(const char* logline, ...);
void debug(const char* logline, ...);
void debug(std::string &str);
void log(std::string &str);
void log(const char* logline, ...);
void error(const char* logline, ...);
void echo (const char* logline, ...);
void echolog(const char* logline, ...);
uint32 tick();

// template class to do convert from
// different types (mainly numerical) 
// to string. It should be just a wrapper
// on C/C++ capabilities of conversation.
// There are two options:
// 1. use stringstream (C++ 03)
// 2. use std::to_string() (C++ 11)
// for more read here:
// http://stackoverflow.com/questions/332111/how-do-i-convert-a-double-into-a-string-in-c
template <class T>
std::string toStr (const T& param) {
  std::stringstream ss;
  ss << param;
  return ss.str();
}

int32 npow(int16 dig, uint16 power);

std::vector<std::string> read_tokens(char *input);
void del_spaces (std::string &str);

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

// this function is used by FEStorage::read_ans_data funciton
uint16 str2dof (std::string dof_name);

char* getCmdOption (char** begin, char** end, const std::string& option);
bool cmdOptionExists (char** begin, char** end, const std::string& option);
std::vector<char*> getCmdManyOptions (char** begin, char** end, const std::string& option); 

struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};

std::string getFileNameFromPath(const std::string filename);

} // namespace nla3d

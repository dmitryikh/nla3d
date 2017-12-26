// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
// before any include we need to define NOMINMAX to avoid redefining `max` with macros
#ifdef WIN32
  #define NOMINMAX
#endif
#include <string>
#include <iostream>
#include <fstream>
#include <strstream>
#include <stdarg.h>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <vector>
#include <set>
#include <list>
#include <algorithm>
#include <stdlib.h>

#define ELPP_STL_LOGGING
#include "easylogging++.h"

#define _USE_MATH_DEFINES
#include <math.h>

#ifdef linux
  #include <string.h>
#endif

// MSVC doent's support C99 standard
#ifdef linux
  #define sprintf_s snprintf
#endif

#ifdef __APPLE__
  #define sprintf_s snprintf
#endif

// TODO:
// for some reasons CHECK_NOTNULL from easylogging++ doesn't work on Apple's clang and g++
#undef CHECK_NOTNULL
#define CHECK_NOTNULL(x) x

// usefull macros to check floats with threshold
#define CHECK_EQTH(a, b, th) CHECK(fabs(a-b) < th)


#define int8 char //-127 to +127
#define uint8 unsigned char //0 to +255
#define int16 short //-32 767 to +32 767
#define uint16 unsigned short //0 to +65 535
#define int32 int //-2 147 483 647 to +2 147 483 647
#define uint32 unsigned int // 0 to +4 294 967 295
#define int64 long long //-9 223 372 036 854 775 807 to +9 223 372 036 854 775 807
#define uint64 unsigned long long //0 to +18 446 744 073 709 551 615


namespace nla3d {

const char SYS_VERSION[] = "1.3";
const char SYS_DATA[] = "16.03.15";

// singleton for configuring the logger
class LogInitializer {
public:
  LogInitializer () {
    el::Configurations conf;
    conf.setToDefault();
    conf.setGlobally(el::ConfigurationType::Format, "%datetime{%H:%m:%s.%g} [%level] %msg");
    conf.setGlobally(el::ConfigurationType::Enabled, "true");
    conf.setGlobally(el::ConfigurationType::ToFile, "true");
    conf.setGlobally(el::ConfigurationType::Filename, "nla3d.log");
    conf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");
    conf.setGlobally(el::ConfigurationType::MillisecondsWidth, "3");
    conf.setGlobally(el::ConfigurationType::PerformanceTracking, "true");
    conf.setGlobally(el::ConfigurationType::MaxLogFileSize, "20971520"); // 20 MB
    conf.setGlobally(el::ConfigurationType::LogFlushThreshold, "10");

    conf.set(el::Level::Debug, el::ConfigurationType::Format, "%datetime{%H:%m:%s.%g} [%level] %func %msg");

    conf.set(el::Level::Error, el::ConfigurationType::Format, "%datetime{%H:%m:%s.%g} [%level] %func %msg");

    // reconfigure all loggers
    el::Loggers::reconfigureAllLoggers(conf);
    el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
  };
};

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

class Timer {
public:
  Timer(bool _start = false) : start_time(0), end_time(0) {
    if (_start) {
      start();
    }
  }

  void start() {
    start_time = clock();
    end_time = start_time;
  }

  double stop() {
    end_time = clock();
    return time();
  }

  double time() {
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

struct MatchPathSeparator {
    bool operator()( char ch ) const {
        return ch == '\\' || ch == '/';
    }
};

std::string getFileNameFromPath(const std::string filename);

enum class ElementType {
  TRUSS3 = 0,
  PLANE41,
  SOLID81,
  TETRA0,
  TETRA1,
  QUADTH,
  SurfaceLINETH,
  TRIANGLE4,
  UNDEFINED
};


const char* const elTypeLabels[] = {
  "TRUSS3",
  "PLANE41",
  "SOLID81",
  "TETRA0",
  "TETRA1",
  "QUADTH",
  "SurfaceLINETH",
  "TRIANGLE4",
  "UNDEFINED"
};

static_assert((int)ElementType::UNDEFINED == sizeof(elTypeLabels)/sizeof(elTypeLabels[0]) - 1,
    "ElementType enumeration and elTypeLabels must have the same number of entries");

} // namespace nla3d

// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <algorithm>
#include <regex>

#include "sys.h"
#include "Node.h"
#include "elements/element.h"
#include "FEStorage.h"
#include "FESolver.h"


using namespace nla3d;


// neutral data structure to keep FE data.
class MeshData {
  public:
    std::vector<uint32> cellNumbers;
    std::vector<std::vector<uint32>> cellNodes;

    // custom element data of different types
    // NOTE: every element in MeshData should have the same set of custom data
    std::map<std::string, std::vector<double>> cellDoubleData;
    std::map<std::string, std::vector<uint32>> cellIntData;
    std::map<std::string, std::vector<std::string>> cellStringData;

    std::vector<uint32> nodesNumbers;
    std::vector<Vec<3>> nodesPos;

    std::vector<loadBC> loadBcs;
    std::vector<fixBC> fixBcs;

    std::vector<Mpc*> mpcs; 

    std::map<std::string, FEComponent> feComps;

    void clear();
    // compress node numbers, we don't compress element numbers because FEStorage will renumber
    // elements anyway
    void compressNumbers();
    // delete nodes which don't belong to any element
    // void deleteOrphanNodes();
    // cell is degenerated when some it's nodes have the same numbers
    std::vector<uint32> getDegeneratedCells();
    std::vector<uint32> getCellsByAttribute(std::string atr_name, uint32 atr_val);

  private:

};

std::string& strim(std::string& str);
// std::move semantic to use with rvalues
std::string&& strim(std::string&& str);
// return first non-space char
char sfirstNotBlank(const std::string& str);

std::string& stolower(std::string& str);
std::string& stoupper(std::string& str);

// split `line` to substrings with `widths` lengths, if `strict` is false then last substrings could
// be missed
std::vector<std::string> ssplit(const std::string& line, const std::vector<int>& widths,
                                  bool strict = true);

// case insensitive string comparison
bool iequals(const std::string& a, const std::string& b);

// get line from the `is`, supported all types of line endings
std::istream& getLine(std::istream& is, std::string& t);


// class to split line into substrings with delimiters.
class Tokenizer {
  public:
    int tokenize(const std::string& line);
    int tokenInt(size_t ind);
    double tokenDouble(size_t ind);
    std::string& getTokenString(size_t ind);

    std::set<char> delimiters;
    std::vector<std::string> tokens;
    std::vector<char> dtokens;
    bool tolower = false;
};


// read Ansys Mechanical APDL *.cdb file. Nodes, Elements, Displacement BC, Force BC and MPC
// (Constraint equations) is supported
bool readCdbFile(std::string filename, MeshData& md);

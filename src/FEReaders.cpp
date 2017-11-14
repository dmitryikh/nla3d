// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "FEReaders.h"
#include <unordered_map>
#include <algorithm>

using std::ifstream;
using std::vector;
using std::string;
using std::map;
using std::cout;
using std::set;
using std::endl;

const string space_chars = " \t";


void MeshData::clear() {
    cellNumbers.clear();
    cellNodes.clear();
    cellDoubleData.clear();
    cellIntData.clear();
    cellStringData.clear();
    nodesNumbers.clear();
    nodesPos.clear();
    loadBcs.clear();
    fixBcs.clear();
    mpcs.clear(); 
    feComps.clear();
}


void MeshData::compressNumbers() {
  // numbers should start from 1
  if (nodesNumbers.size() == *std::max_element(nodesNumbers.begin(), nodesNumbers.end())) { 
    // nothing to do
    return;
  }
  std::map<uint32, uint32> old2new;
  uint32 nextNumber = 1;
  for (uint32 i = 0; i < nodesNumbers.size(); i++) {
    old2new[nodesNumbers[i]] = nextNumber;
    nodesNumbers[i] = nextNumber;
    nextNumber++;
  }
  assert(nextNumber - 1 == nodesNumbers.size());

  // go through elements and change element's node numbers
  for (auto& e : cellNodes) {
    for (uint16 i = 0; i < e.size(); i++) {
      e[i] = old2new[e[i]];
    }
  }

  // go through loadBCs
  for (auto& bc : loadBcs) {
    bc.node = old2new[bc.node];
  }

  // go through fixBCs
  for (auto& bc : fixBcs) {
    bc.node = old2new[bc.node];
  }

  // go through Mpcs
  for (auto& mpc : mpcs) {
    for (auto& term : mpc->eq) {
      term.node = old2new[term.node];
    }
  }

  // go through FE components
  for (auto& pair : feComps) {
    auto& comp = pair.second;
    if (comp.type == FEComponent::NODES) {
      for (auto& n : comp.list) {
        n = old2new[n];
      } 
    }
  }
}


std::vector<uint32> MeshData::getDegeneratedCells() {
  std::vector<uint32> res;
  for (size_t i = 0; i < cellNumbers.size(); i++) {
    std::set<uint32> s(cellNodes[i].begin(), cellNodes[i].end());
    if (s.size() != cellNodes[i].size()) {
      res.push_back(i);
    }
  }
  return res;
}


std::vector<uint32> MeshData::getCellsByAttribute(std::string atr_name, uint32 atr_val) {
  std::vector<uint32> res;
  for (size_t i = 0; i < cellNumbers.size(); i++) {
    if (cellIntData[atr_name][i] == atr_val) {
      res.push_back(i);
    }
  }
  return res;
}


string&& strim(string&& str) {
  size_t st = str.find_first_not_of(space_chars, 0);
  str.erase(0, st);
  size_t en = str.find_last_not_of(space_chars);
  str.erase(en+1);
  return std::move(str);
}


string& strim(string& str) {
  size_t st = str.find_first_not_of(space_chars, 0);
  str.erase(0, st);
  size_t en = str.find_last_not_of(space_chars);
  str.erase(en+1);
  return str;
}


string& stolower(string& str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  return str;
}


string& stoupper(string& str) {
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
  return str;
}


char sfirstNotBlank(const string& str) {
  size_t st = str.find_first_not_of(space_chars, 0);
  return str[st];
}


std::vector<std::string> ssplit(const std::string& line, const std::vector<int>& widths, bool strict) {
  int ind = 0;
  vector<string> vv;
  for (int w : widths) {
    if (line.length() - ind < w) {
      // no room for another field
      if (strict) {
        LOG(FATAL) << "ssplit: incomplete line: \"" << line << "\"";
      } else {
        if (line.length() - ind != 0) {
          LOG(FATAL) << "ssplit: incomplete field: \"" << line << "\"";
        }
        return std::move(vv);
      }
    }
    vv.push_back(line.substr(ind, w));
    ind += w;
  }
  return vv;
}

 
bool iequals(const string& a, const string& b) {
    unsigned int sz = a.size();
    if (b.size() != sz)
        return false;
    for (unsigned int i = 0; i < sz; ++i)
        if (tolower(a[i]) != tolower(b[i]))
            return false;
    return true;
}

// taken from
// http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf/6089413#6089413 
std::istream& getLine(std::istream& is, std::string& t) {
  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  std::istream::sentry se(is, true);
  std::streambuf* sb = is.rdbuf();

  for (;;) {
    int c = sb->sbumpc();
    switch (c) {
      case '\n':
        return is;
      case '\r':
        if (sb->sgetc() == '\n') {
          sb->sbumpc();
        }
        return is;
      case EOF:
        // Also handle the case when the last line has no line ending
        if (t.empty()) {
          is.setstate(std::ios::eofbit);
        }
        return is;
      default:
        t += (char)c;
    }
  }
}


int Tokenizer::tokenize(const string& line) {
  tokens.clear();
  dtokens.clear();

  string _delimiters(delimiters.begin(), delimiters.end());

  int found = 0;
  size_t st = 0;

  while (1) {
    size_t ind = line.find_first_of(_delimiters, st);
    found++;
    if (ind == string::npos) {
      tokens.push_back(strim(line.substr(st, string::npos)));
      break;
    } else {
      tokens.push_back(strim(line.substr(st, ind - st)));
      // if we want to store which exactly delimeter was used
      dtokens.push_back(line[ind]);
    }
    st = ind+1;
  }

  if (tolower) {
    for (int i = 0; i < tokens.size(); i++) {
      stolower(tokens[i]);
    }
  }
  return found;
}


int Tokenizer::tokenInt(size_t ind) {
  return std::stoi(tokens[ind]);
}


double Tokenizer::tokenDouble(size_t ind) {
  return std::stod(tokens[ind]);
}


bool readCdbFile(std::string filename, MeshData& md) {
  uint32 n_number, en;
  ifstream file(filename);
  if (!file) {
    LOG(WARNING) << "Can't open cdb file " << filename;
    return false;
  }
  md.clear();
  // init element addition data which is specific to cdb file format
  md.cellIntData.insert(std::make_pair("MAT", std::vector<uint32>()));
  md.cellIntData.insert(std::make_pair("REAL", std::vector<uint32>()));
  md.cellIntData.insert(std::make_pair("SECT", std::vector<uint32>()));
  md.cellIntData.insert(std::make_pair("CS", std::vector<uint32>()));
  md.cellIntData.insert(std::make_pair("SHAPE", std::vector<uint32>()));
  md.cellIntData.insert(std::make_pair("TYPE", std::vector<uint32>()));

  std::unordered_map<std::string, double> apdlParameters;

  // current element type number. We keep the current element type number while proceed the apld
  // file. Currently, this info is not used, but in the past it was used to determine element type
  // while parsing EBLOCK blocks. 
  uint32 type = 0;

  Tokenizer t;
  // do not convert tokens to lower symbols
  t.tolower = false;
  t.delimiters.insert(',');

  string line;
  while (!getLine(file, line).eof()) {
    // split line with delimiters ','
    t.tokenize(line);
    if (t.tokens[0] == "") continue;
    if (t.tokens[0][0] == '/') continue;
    if (t.tokens[0][0] == '!') continue;

    if (iequals(t.tokens[0], "NBLOCK")) {
      // example from cdb Ansys APDL 12.1:
      //NBLOCK,6,SOLID,     9355,     9355
      //(3i9,6e20.13)
      //        1        0        0 7.0785325971794E+01 6.5691449317818E+01-3.6714639015390E+01
      //
      // from Ansys WB 17.1:
      // nblock,3
      // (1i9,3e20.9e3)
      //         1    2.000000000E+000    6.000000000E+000    0.000000000E+000
      //         2    0.000000000E+000    6.000000000E+000    0.000000000E+000
      //
      // from Ansys WB APDL 17.1: 
      //NBLOCK,6,SOLID,       341,       341
      //(3i9,6e21.13e3)
      //        1        0        0 0.0000000000000E+000 4.0000000000000E+000
      //        2        0        0 0.0000000000000E+000 6.0000000000000E+000
      // parse format string
      getLine(file, line);
      stolower(strim(line));
      std::regex re(R"(\((\d+)i(\d+),(\d+)e(\d+)\.[0-9e]+\))");
      std::smatch match;
      int int_num;
      int int_field;
      int float_num;
      int float_field;
      if (std::regex_match(line, match, re)) {
        int_num = std::stoi(match[1]);
        int_field = std::stoi(match[2]);
        float_num = std::stoi(match[3]);
        float_field = std::stoi(match[4]);
      } else {
        LOG(FATAL) << "Don't understand nblock format string \"" << line << "\"";
      }

      // prepare array of field widths
      std::vector<int> widths;
      for (int i = 0; i < int_num; i++) {
        widths.push_back(int_field);
      }
      for (int i = 0; i < float_num; i++) {
        widths.push_back(float_field);
      }

      // read nodes info line by line
      while(!getLine(file, line).eof()) {
        // check for the end of nodal table
        char ch = sfirstNotBlank(line);
        if (ch == '-' || ch == 'N') break;

        auto v = ssplit(line, widths, false);
        // column no. 1 - node number
        uint32 number = std::stoi(v[0]);
        // Vec is initialized with zeros by default
        Vec<3> pos;
        // cdb file could skip last components of positions if they are equal to zero
        if (v.size() > int_num) pos[0] = std::stod(v[int_num]);
        if (v.size() > int_num + 1) pos[1] = std::stod(v[int_num + 1]);
        if (v.size() > int_num + 2) pos[2] = std::stod(v[int_num + 2]);
        md.nodesNumbers.push_back(number);
        md.nodesPos.push_back(pos);
      }
    }//NBLOCK
    else if (iequals(t.tokens[0], "EBLOCK")) {  
      // Example of record:
      //EBLOCK,19,SOLID,      7024,      7024
      //(19i9)
      // ...
      // Info from Ansys 17 help pages:
      /*The format of the element "block" is as follows for the SOLID format:
      Field 1 - The material number.
      Field 2 - The element type number.
      Field 3 - The real constant number.
      Field 4 - The section ID attribute (beam section) number.
      Field 5 - The element coordinate system number.
      Field 6 - The birth/death flag.
      Field 7 - The solid model reference number.
      Field 8 - The element shape flag.
      Field 9 - The number of nodes defining this element if Solkey = SOLID; otherwise, Field 9 = 0.
      Field 10 - Not used.
      Field 11 - The element number.
      Fields 12-19 - The node numbers. The next line will have the additional node numbers if there
      are more than eight.

      The format without the SOLID keyword is:
      Field 1 - The element number.
      Field 2 - The type of section ID.
      Field 3 - The real constant number.
      Field 4 - The material number.
      Field 5 - The element coordinate system number.
      Fields 6-15 - The node numbers. The next line will have the additional node numbers if there
      are more than ten.
      The final line of the block will be a -1 in field 1. */

      //TODO: It seems that getline keeps windows line ending
      bool isSolid = iequals(t.tokens[2], "SOLID");
      // parse format string
      getLine(file, line);
      stolower(strim(line));
      std::regex re(R"(\((\d+)i(\d+)\))");
      std::smatch match;
      int int_num;
      int int_field;
      if (std::regex_match(line, match, re)) {
        int_num = std::stoi(match[1]);
        int_field = std::stoi(match[2]);
      } else {
        LOG(FATAL) << "Don't understand eblock format string \"" << line << "\"";
      }

      // prepare array of field widths
      std::vector<int> widths;
      for (int i = 0; i < int_num; i++) {
        widths.push_back(int_field);
      }
      // read element data line by line
      while(!getLine(file, line).eof()) {
        // check for the end of element table
        char ch = sfirstNotBlank(line);
        if (ch == '-' || ch == 'N') break;

        auto v = ssplit(line, widths, false);

        // initialize values that we read from element table
        uint32 MAT = 0;
        uint32 ETYPE = 0;
        uint32 REAL = 0;
        uint32 SECT = 0;
        uint32 CS = 0;
        uint32 SHAPE = 0;
        uint32 nNodes = 0;
        uint32 elNumber = 0;
        size_t st = 0;
        if (isSolid) {
          MAT = static_cast<uint32>(std::stoi(v[0]));
          ETYPE = static_cast<uint32>(std::stoi(v[1]));
          REAL = static_cast<uint32>(std::stoi(v[2]));
          SECT = static_cast<uint32>(std::stoi(v[3]));
          CS = static_cast<uint32>(std::stoi(v[4]));
          SHAPE = static_cast<uint32>(std::stoi(v[7]));
          nNodes = static_cast<uint32>(std::stoi(v[8]));
          elNumber = static_cast<uint32>(std::stoi(v[10]));
          st = 11;
        } else {
          elNumber = static_cast<uint32>(std::stoi(v[0]));
          // based on my experiments on *.apdl files this field is about etype, but not about secnum
          ETYPE = static_cast<uint32>(std::stoi(v[1]));
          REAL = static_cast<uint32>(std::stoi(v[2]));
          MAT = static_cast<uint32>(std::stoi(v[3]));
          CS = static_cast<uint32>(std::stoi(v[4]));

          // in this case we need to determine number of nodes by ourselfs
          nNodes = v.size() - 5;
          st = 5;
        }
        vector<uint32> enodes;
        
        if (v.size() != st + nNodes) {
          LOG(FATAL) << "Not enought fields to read element nodes";
        }
        // parse element nodes
        for (int i = 0; i < nNodes; i++) {
          enodes.push_back(static_cast<uint32>(std::stoi(v[st + i])));
        }

        md.cellNumbers.push_back(elNumber);
        md.cellNodes.push_back(enodes);
        md.cellIntData["TYPE"].push_back(ETYPE);
        md.cellIntData["MAT"].push_back(MAT);
        md.cellIntData["REAL"].push_back(REAL);
        md.cellIntData["SECT"].push_back(SECT);
        md.cellIntData["CS"].push_back(CS);
        md.cellIntData["SHAPE"].push_back(SHAPE);
      } // end of element table
    }//EBLOCK
    else if (iequals(t.tokens[0], "D")) {
      // Fixed dof record
      fixBC bnd;
      bool isNumeric = false;
      try {
        bnd.node = t.tokenInt(1);
        bnd.value = t.tokenDouble(3);
        isNumeric = true;
      } catch (const std::invalid_argument& ia) {
        isNumeric = false;
      }
      if (isNumeric) {
        // node and value field are number, we can deal with in. In other cases these can be 'all',
        bnd.node_dof = Dof::label2dofType(t.tokens[2]);
        md.fixBcs.push_back(bnd);
      }
    }//D
    else if (iequals(t.tokens[0], "CE")) {
      // CE - is a command to declare MPC equation
      // How MPC looks like this in cdb file:
      //CE,R5.0,DEFI,       2,       1,  0.00000000    
      //CE,R5.0,NODE,      1700,UX  ,  1.00000000    ,      1700,UZ  ,  1.00000000  
      // TODO: Mpc is dynamicaly allocated, but MeshData won't free this memory (this should be done
      // in FEStorage). This is potential memory leak
      Mpc* mpc = new Mpc();
      mpc->b = t.tokenDouble(5);
      int n_terms = t.tokenInt(3);
      // read MPC terms. They are stored by 2 in a row
      while (n_terms > 0) {
        getLine(file, line);
        t.tokenize(line);
        uint16 place = 6;
        for (int i = 0; i < std::max(n_terms, 2); i++) {
          uint32 node = t.tokenInt(3 + 3 * i + 0);
          Dof::dofType dof = Dof::label2dofType(t.tokens[3 + 3 * i + 1]);
          double coef = t.tokenDouble(3 + 3 * i + 2);
          mpc->eq.push_back(MpcTerm(node,dof,coef));
          n_terms--;
        }
      }
      md.mpcs.push_back(mpc);
    }//CE (MPC)
    else if (iequals(t.tokens[0], "CMBLOCK")) {
      // Example of CMBLOCK command:
      //CMBLOCK,BOTTOM_SIDE,NODE,      17  ! users node component definition
      //(8i10)
      //      5037     -5341      6330     -6352      6355      6357      6433     -6456
      //      6459      6470     -6473      6537     -6556      6566     -6569      6633
      //     -6652
      if (t.tokens[1][0] != '_') {
        // we work only with component names started not with underscore '_'
        FEComponent comp;
        comp.name = t.tokens[1];
        // NODE or ELEMENT component?
        comp.type = FEComponent::typeFromString(stoupper(t.tokens[2]));

        size_t numRanges = t.tokenInt(3);
        getLine(file, line);
        // we need to take a format of columns "(8i10)"
        // read how many records in a single row
        stolower(strim(line));
        std::regex re(R"(\((\d+)i(\d+)\))");
        std::smatch match;
        int int_num;
        int int_field;
        if (std::regex_match(line, match, re)) {
          int_num = std::stoi(match[1]);
          int_field = std::stoi(match[2]);
        } else {
          LOG(FATAL) << "Don't understand CMBLOCK format string \"" << line << "\"";
        }
        // prepare array of field widths
        std::vector<int> widths;
        for (int i = 0; i < int_num; i++) {
          widths.push_back(int_field);
        }

        uint16 in_row = 0;
        uint16 all = 0;

        // read all ranges and calculate overall number of entities in component list
        vector<int32> rangesVec;
        rangesVec.reserve(numRanges);
        getLine(file, line);
        auto v = ssplit(line, widths, false);
        size_t numEntity = 0;
        while (all < numRanges) {
           if (in_row == int_num) {
              getLine(file, line);
              v = ssplit(line, widths, false);
              in_row = 0;
           }
           rangesVec.push_back(std::stoi(v[in_row]));
           if (rangesVec[all] > 0) {
              numEntity++;
           } else {
              // if next index is negative that means a range with step 1
              // for example: 4 -9 --> 4 5 6 7 8 9
              assert(all > 0);
              assert(rangesVec[all-1] > 0);
              numEntity += -rangesVec[all] - rangesVec[all-1];
           }
           in_row ++;
           all ++;
        }
        // fill list of entities numbers
        comp.list.reserve(numEntity);
        all = 0;
        while (all < numRanges) {
          if (rangesVec[all] > 0) {
            comp.list.push_back(rangesVec[all]);
          } else {
            for (uint32 i = rangesVec[all-1] + 1; i < static_cast<uint32> (-rangesVec[all]+1); i++) {
              // TODO: need assert for overflow
              comp.list.push_back(i);
            }
          }
          all++;
        }
        md.feComps.insert(std::make_pair(comp.name, comp));
      }
    } //CMBLOCK
    else if (iequals(t.tokens[0], "ET") || iequals(t.tokens[0], "ETYPE")) {
      // Keep a track on last ET command.. We need this info to know element type.
      // If we can convert ET argument into int, then we try to find the APDL parameter in
      // apdlParameters dictionary. If it's not there - adpParameters will return new zero value.
      // That means, that if parameter is not declared(somehow..) we just put type = 0
      // NOTE: in current implementation type value is not used, but it could be used again in the
      // future.
      try {
        type = t.tokenInt(1);
      } catch (const std::invalid_argument& ia) {
        type = static_cast<uint32>(apdlParameters[t.tokens[1]]);
      }
    } // ETYPE
    else if (iequals(t.tokens[0], "TYPE")) {
      // Keep the track on last TYPE command.. We need this info to know element type.
      // NOTE: in current implementation type value is not used, but it could be used again in the
      // future.
      try {
        type = t.tokenInt(1);
      } catch (const std::invalid_argument& ia) {
        type = 0;
      }
    } // ETYPE
    else if (iequals(t.tokens[0], "F")) {
      // Dof force record
      //f,42,heat,800000.
      loadBC bnd;
      bool isNumeric = false;
      try {
        bnd.node = t.tokenInt(1);
        bnd.value = t.tokenDouble(3);
        isNumeric = true;
      } catch (const std::invalid_argument& ia) {
        isNumeric = false;
      }
      if (isNumeric) {
        // node and value field are number, we can deal with in. In other cases these can be 'all',
        // parameter or reference to table '%my_apdl_table%. We can't work with these complex cases.
        if (iequals(t.tokens[2], "HEAT")) {
          bnd.node_dof = Dof::TEMP;
        } else if (iequals(t.tokens[2], "FX")) {
          bnd.node_dof = Dof::UX;
        } else if (iequals(t.tokens[2], "FY")) {
          bnd.node_dof = Dof::UY;
        } else if (iequals(t.tokens[2], "FZ")) {
          bnd.node_dof = Dof::UZ;
        } else {
          bnd.node_dof = Dof::label2dofType(t.tokens[2]);
        }
        md.loadBcs.push_back(bnd);
      }
    }//F
    else if (iequals(t.tokens[0], "*set")) {
      // This is declaration of APDL parameters. We track all APDL numeric variables which is
      // defined in APDL file.
      // NOTE: currently this information is not used, but it could be useful later
      bool isNumeric = false;
      double value;
      std::string parName;
      // try to read parameter's name and it's value as double
      // if fail - skip this expression
      try {
        parName = t.tokens[1];
        value = t.tokenDouble(2);
        isNumeric = true;
      } catch (const std::invalid_argument& ia) {
        isNumeric = false;
      }
      if(isNumeric) {
        apdlParameters[parName] = value;
      }
    }
  }
  file.close();
  return true;
}

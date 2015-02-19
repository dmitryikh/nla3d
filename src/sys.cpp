#include "sys.h"
#include "FE_Storage.h"

static Log_opts log_opts;


void warning(const char* logline, ...)
{
	WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList); 
	va_end(argList);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << "Warning: " << buffer <<endl;
		ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout << "Warning: " << buffer <<endl;
	file << "Warning: " << buffer <<endl;
	file.close();
	ReleaseMutex( log_opts.output_lock );
}

void debug(const char* logline, ...)
{
	if (!debug_mode) return;
	WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	va_end(argList);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << "DEBUG: " << buffer <<endl;
		ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout  << "DEBUG: "<< buffer <<endl;
	file  << "DEBUG: "<< buffer <<endl;
	file.close();
	ReleaseMutex( log_opts.output_lock );
}

void debug (string &str)
{
	WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout <<  str <<endl;
		ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout  << "DEBUG: "<< str <<endl;
	file  << "DEBUG: "<< str <<endl;
	file.close();
	ReleaseMutex( log_opts.output_lock );

}

void log(string &str)
{
	WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << str <<endl;
		ReleaseMutex( log_opts.output_lock );
		return;
	}
	file << str <<endl;
	file.close();
	ReleaseMutex( log_opts.output_lock );
}

void log(const char* logline, ...)
{
	WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	va_end(argList);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << buffer <<endl;
		ReleaseMutex( log_opts.output_lock );
		return;
	}
	file << buffer <<endl;
	file.close();
	ReleaseMutex( log_opts.output_lock );
}

void error(const char* logline, ...)
{
	WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	va_end(argList);
	if (file)
	{
		file << "Error: " << buffer <<endl;
		file.close();
	}
	ReleaseMutex( log_opts.output_lock );
	exit(0);
}

void echo(const char* logline, ...)
{
	WaitForSingleObject( log_opts.output_lock, INFINITE );
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	va_end(argList);
	cout << buffer <<endl;
	ReleaseMutex( log_opts.output_lock );
}

void echolog(const char* logline, ...)
{
	WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	va_end(argList);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << buffer <<endl;
		ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout << buffer <<endl;
	file << buffer <<endl;
	file.close();
	ReleaseMutex( log_opts.output_lock );
}

uint32 tick()
{
	return (uint32) clock();
}

string IntToStr (uint32 dig)
{
	stringstream ss;
	ss << dig;
	return ss.str();
}

int32 npow(int16 dig, uint16 power)
{
	int32 res=1; //TODO: if too big number?
	for (uint16 i=0; i < power; i++)
		res = res * dig;
	return res;
}

vector<string> read_tokens(char *input)
{
	vector<string> vec;
	string tmp("");
	char delimeters[]="(),";
	char *p=input;
	char *start=p;
	while (*p)
	{
		if (*p=='!') break;
		bool isfind=false;
		char* delp=delimeters;
		while (*delp)
		{
			if (*delp==*p)
			{
				isfind=true;
				break;
			}
			delp++;
		}
		if (isfind)
		{
			tmp.assign(start, (int16) (p-start));
			del_spaces(tmp);
			vec.push_back(tmp);
			vec.push_back(string(delp,1));
			p++;
			start=p;
			continue;
		}
		p++;
	}
	tmp.assign(start, (int16) (p-start));
	del_spaces(tmp);
	vec.push_back(tmp);
	return vec;
}

void del_spaces (string &str)
{
	uint16 start = 0;
	if (str.length()==0) return;
	while (str[start]==' ' ) 
	{
		start++;
		if (start == str.length()) 
		{
			str=" ";
			return;
		}
	}
	uint16 end = str.length()-1;
	while (str[end] ==' ') end--;
	str=string(str,start, end-start+1);
}

// read Ansys Mechanical APDL *.cdb file. Nodes, Elements, Displacement BC and MPC (Constraint equations) is supported
// read_ans_data repcales storage's mesh.
bool read_ans_data(const char *filename, FE_Storage_Interface *storage)
{
	uint32 n_number, en;
	ifstream file(filename);
	if (!file)
	{
		warning("read_ans_data: Can't open input file `%s`",filename);
		return false;
	}
	storage->clearMesh();
	char buf[1024]="";
	while (file.getline(buf, 1024))
	{
		vector<string> vec = read_tokens(buf);
		if (vec[0].compare("NBLOCK") == 0)
		{
    //NBLOCK,6,SOLID,     9355,     9355
    //(3i9,6e20.13)
    //        1        0        0 7.0785325971794E+01 6.5691449317818E+01-3.6714639015390E+01
			uint32 max_n_number = atoi(vec[6].c_str());
			n_number= atoi(vec[6].c_str());
      if (max_n_number != n_number) {
        warning("read_ans_data: NBLOCK: maximum node number is %d, but number of nodes is %s. Note that nla3d needs compressed numbering for nodes and elements", max_n_number, n_number );
      }
			storage->nodes_reassign(n_number);
			file.getline(buf, 1024);
      string buf_str(buf);
      // we need to take a format of columns "3i9"
      uint16 start = buf_str.find("i")+1;
      uint16 stop = buf_str.find(",");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str());

      //TODO: as we need numbering from 1/0 to n_number, here we can check that number of nodes and theirs id are in a row
			for (uint32 i=1; i<=n_number; i++)
			{
				file.getline(buf, 1024);
				uint16 len=strlen(buf);
				for (uint16 j=0; j<3;j++)
					if (len>=3*frmt+20*(j+1))
            //note that last column in NBLOCK table could be avoided if Z=0.0
            //but storage->nodes_reassign(n_number) initialize the node table with (0,0,0)
						storage->getNode(i).pos[j] = atof(string((char*) (buf+3*frmt+20*j),20).c_str());
			}
		}//NBLOCK
		else if (vec[0].find("EBLOCK")!=vec[0].npos)
		{  
      //EBLOCK,19,SOLID,      7024,      7024
      //(19i9)
			en=atoi(vec[6].c_str());
			storage->elements_reassign(en);
			file.getline(buf, 1024);
      // we need to take a format of columns "3i9"
      // in Ansys 12 here is 8 symbols per number (19i8), but in ansys 15 (19i9) is used. 
      string buf_str(buf);
      uint16 start = buf_str.find("i")+1;
      uint16 stop = buf_str.find(")");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str()); 
			for (uint32 i=1; i<=en; i++)
			{
				file.getline(buf, 1024);
				uint16 len=strlen(buf);
        //TODO: check that n_nodes and provided number of nodes are the same
        if (len != 11*frmt+frmt*Element::n_nodes())
          warning("read_ans_data: in EBLOCK for element %d the number of nodes provided is not equal to %d", i, Element::n_nodes());
				for (uint16 j=0; j<Element::n_nodes();j++)
					if (len>=11*frmt+frmt*(j+1))
						storage->getElement(i).node_num(j) = atoi(string((char*) (buf+11*frmt+frmt*j),frmt).c_str());
			}
		}//EBLOCK
		else if (vec[0].find('D')!=vec[0].npos && vec[0].length()==1)
		{
				BC_dof_constraint bnd;
				bnd.node = atoi(vec[2].c_str());
				bnd.value = atof(vec[6].c_str());
        bnd.node_dof = str2dof(vec[4]);
        storage->add_bounds(bnd);
		}//D
    else if (vec[0].compare("CE") == 0)
    {
      //How MPC looks like this in cdb file:
      //CE,R5.0,DEFI,       2,       1,  0.00000000    
      //CE,R5.0,NODE,      1700,UX  ,  1.00000000    ,      1700,UZ  ,  1.00000000  
      BC_MPC mpc;
      mpc.b = atoi(vec[10].c_str()); //rhs of MPC equation
      uint16 n_terms = atoi(vec[6].c_str()); //number of terms in equation
      //debug("MPC link: %d terms, b = %f", n_terms, mpc.b);
      while (n_terms > 0)
      {
        file.getline(buf, 1024);
        vector<string> vec = read_tokens(buf);
        uint16 place = 6;
        for (int i=0; i < max(n_terms, 2); i++) 
        {
          uint32 node = atoi(vec[place].c_str());
          uint16 dof = str2dof(vec[place+2]);
          double coef = atof(vec[place+4].c_str());
          //debug("%d term: node %d, dof %d, coef = %f", i, node, dof, coef);
          mpc.eq.push_back(MPC_token(node,dof,coef));
          place += 6;
          n_terms--;
        }
      }
			storage->add_bounds(mpc);
    }//CE (MPC)
    //TODO: add FX FY FZ
	}
	file.close();
	storage->setStatus(ST_LOADED);
	return true;
}

//TODO: use Doftype
uint16 str2dof (string dof_name)
{
  uint16 dof;
  if (!dof_name.compare("UX"))
  {
    dof = 0;
  }
  else if (!dof_name.compare("UY"))
  {
    dof = 1;
  }
  else if (!dof_name.compare("UZ"))
  {
    dof = 2;
  }
  else
  {
    warning("str2dof: unknown dof key: %s", dof_name);
    dof = 0; //TODO: what we need to return in this case??
  }
  return dof;
}


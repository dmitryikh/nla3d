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
		cout << "Warning: Can't open log file: " << log_file_name << endl;
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


bool read_ans_data(const char *filename, FE_Storage_Interface *storage)
{
	//функция загружает КЭ из Ansys
	uint32 nn, en;
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
		if (vec[0].find("NBLOCK")!=vec[0].npos)
		{
			nn = atoi(vec[6].c_str());
			storage->nodes_reassign(nn);
			file.getline(buf, 1024);
			for (uint32 i=1; i<=nn; i++)
			{
				file.getline(buf, 1024);
				uint16 len=strlen(buf);
				for (uint16 j=0; j<3;j++)
					if (len>=3*8+20*(j+1))
						storage->getNode(i).pos[j] = atof(string((char*) (buf+3*8+20*j),20).c_str());
			}
		}//NBLOCK
		else if (vec[0].find("EBLOCK")!=vec[0].npos)
		{  
			en=atoi(vec[6].c_str());
			storage->elements_reassign(en);
			file.getline(buf, 1024);
			for (uint32 i=1; i<=en; i++)
			{
				file.getline(buf, 1024);
				uint16 len=strlen(buf);
				for (uint16 j=0; j<Element::n_nodes();j++)
					if (len>=11*8+8*(j+1))
						storage->getElement(i).node_num(j) = atoi(string((char*) (buf+11*8+8*j),8).c_str());
			}
		}//EBLOCK
		else if (vec[0].find('D')!=vec[0].npos && vec[0].length()==1)
		{
				BC_dof_constraint bnd;
				bnd.node = atoi(vec[2].c_str());
				bnd.node_dof = 0;
				bnd.value = atof(vec[6].c_str());
				if (!vec[4].compare("UX"))
				{
					bnd.node_dof = 0;
					storage->add_bounds(bnd);
				}
				else if (!vec[4].compare("UY"))
				{
					bnd.node_dof = 1;
					storage->add_bounds(bnd);
				}
				else if (!vec[4].compare("UZ"))
				{
					bnd.node_dof = 2;
					storage->add_bounds(bnd);
				}
				else
				{
					warning("read_ans_data: unknown D key");
				}
				//TODO: add FX FY FZ
		}//D
	}
	file.close();
	storage->setStatus(ST_LOADED);
	return true;
}

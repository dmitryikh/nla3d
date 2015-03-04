#include "sys.h"
#include "FE_Storage.h"

static Log_opts log_opts;


void warning(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList); 
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << "Warning: " << buffer <<endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout << "Warning: " << buffer <<endl;
	file << "Warning: " << buffer <<endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void debug(const char* logline, ...)
{
	if (!debug_mode) return;
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
//	vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << "DEBUG: " << buffer <<endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout  << "DEBUG: "<< buffer <<endl;
	file  << "DEBUG: "<< buffer <<endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void debug (string &str)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout <<  str <<endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout  << "DEBUG: "<< str <<endl;
	file  << "DEBUG: "<< str <<endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );

}

void log(string &str)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << str <<endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	file << str <<endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void log(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << buffer <<endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	file << buffer <<endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void error(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);


	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << "Error: " << buffer <<endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout << "Error: " << buffer <<endl;
	file << "Error: " << buffer <<endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
	exit(1);
}

void echo(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	cout << buffer <<endl;
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void echolog(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	ofstream file(log_file_name,ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	if (!file)
	{
		cout << "Warning: Can't open log file: " << log_file_name <<endl;
		cout << buffer <<endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	cout << buffer <<endl;
	file << buffer <<endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
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
    warning("str2dof: unknown dof key: %s", dof_name.c_str());
    dof = 0; //TODO: what we need to return in this case??
  }
  return dof;
}

char* getCmdOption(char ** begin, char ** end, const std::string & option) {
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return *itr;
    }
    return 0;
}

vector<char*> getCmdManyOptions(char ** begin, char ** end, const std::string & option) {
    vector<char*> _vec;
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
      while(*itr[0] != '-') {
        _vec.push_back(*itr);
        if (++itr == end) {
          break;
        }
      }
    }
    return _vec;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}


string getFileNameFromPath(const string filename) {
    std::string::const_reverse_iterator pivot = 
          std::find( filename.rbegin(), filename.rend(), '.' );
    return pivot == filename.rend()
        ? filename
        : std::string( filename.begin(), pivot.base() - 1 );
}

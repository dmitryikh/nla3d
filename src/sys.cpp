// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "sys.h"
#include "FEStorage.h"

namespace nla3d {

static Log_opts log_opts;

void warning(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	std::ofstream file(log_file_name,std::ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList); 
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	if (!file)
	{
		std::cout << "Warning: Can't open log file: " << log_file_name <<std::endl;
		std::cout << "Warning: " << buffer <<std::endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	std::cout << "Warning: " << buffer <<std::endl;
	file << "Warning: " << buffer <<std::endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void debug(const char* logline, ...)
{
	if (!debug_mode) return;
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	std::ofstream file(log_file_name,std::ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
//	vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	if (!file)
	{
		std::cout << "Warning: Can't open log file: " << log_file_name <<std::endl;
		std::cout << "DEBUG: " << buffer <<std::endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	std::cout  << "DEBUG: "<< buffer <<std::endl;
	file  << "DEBUG: "<< buffer <<std::endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void debug (std::string &str)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	std::ofstream file(log_file_name,std::ios::app);
	if (!file)
	{
		std::cout << "Warning: Can't open log file: " << log_file_name <<std::endl;
		std::cout <<  str <<std::endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	std::cout  << "DEBUG: "<< str <<std::endl;
	file  << "DEBUG: "<< str <<std::endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );

}

void log(std::string &str)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	std::ofstream file(log_file_name,std::ios::app);
	if (!file)
	{
		std::cout << "Warning: Can't open log file: " << log_file_name <<std::endl;
		std::cout << str <<std::endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	file << str <<std::endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void log(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	std::ofstream file(log_file_name,std::ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	if (!file)
	{
		std::cout << "Warning: Can't open log file: " << log_file_name <<std::endl;
		std::cout << buffer <<std::endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	file << buffer <<std::endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void error(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	std::ofstream file(log_file_name,std::ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);


	if (!file)
	{
		std::cout << "Warning: Can't open log file: " << log_file_name <<std::endl;
		std::cout << "Error: " << buffer <<std::endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	std::cout << "Error: " << buffer <<std::endl;
	file << "Error: " << buffer <<std::endl;
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
	std::cout << buffer <<std::endl;
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

void echolog(const char* logline, ...)
{
  //TODO: mutex
	//WaitForSingleObject( log_opts.output_lock, INFINITE );
	std::ofstream file(log_file_name,std::ios::app);
	va_list argList;
	char buffer[1024];
	va_start(argList, logline);
	//vsnprintf_s(buffer, 1024,_TRUNCATE, logline, argList);
	vsnprintf(buffer, 1024, logline, argList); 
	va_end(argList);
	if (!file)
	{
		std::cout << "Warning: Can't open log file: " << log_file_name <<std::endl;
		std::cout << buffer <<std::endl;
    //TODO: mutex
    //ReleaseMutex( log_opts.output_lock );
		return;
	}
	std::cout << buffer <<std::endl;
	file << buffer <<std::endl;
	file.close();
  //TODO: mutex
	//ReleaseMutex( log_opts.output_lock );
}

uint32 tick()
{
	return (uint32) clock();
}

int32 npow(int16 dig, uint16 power)
{
	int32 res=1; //TODO: if too big number?
	for (uint16 i=0; i < power; i++)
		res = res * dig;
	return res;
}

std::vector<std::string> read_tokens(char *input)
{
	std::vector<std::string> vec;
	std::string tmp("");
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
			vec.push_back(std::string(delp,1));
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

void del_spaces (std::string &str)
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
	str=std::string(str,start, end-start+1);
}


//TODO: use Doftype
uint16 str2dof (std::string dof_name)
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

std::vector<char*> getCmdManyOptions(char ** begin, char ** end, const std::string & option) {
    std::vector<char*> _vec;
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
      //TODO: what if len(itr[0]) == 1 ?? momry corrupted?
      while(!(itr[0][0] == '-' && (itr[0][1] < '0' || itr[0][1] > '9'))) {
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

//TODO: this functions only truncate file extension. 
//But it was intended to leave only a file name (delete path and extension)
std::string getFileNameFromPath(const std::string filename) {
    std::string::const_reverse_iterator pivot = 
          std::find( filename.rbegin(), filename.rend(), '.' );
    return pivot == filename.rend()
        ? filename
        : std::string( filename.begin(), pivot.base() - 1 );
}

} // namespace nla3d

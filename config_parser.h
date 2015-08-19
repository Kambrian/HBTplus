#ifndef CONFIG_PARSER_H_INCLUDED
#define CONFIG_PARSER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
using namespace std;

#define NumberOfConfigEntries 5

class Parameter_t
{
public:
  string SnapshotPath;
  string SnapshotFileBase;
  string HaloPath;
  string SubhaloPath;
  int MinNumPartOfSub;
  Parameter_t(): MinNumPartOfSub(10)
  {
  }
  void ParseConfigFile(char * param_file);
  void SetParameterValue(const string &line);
};

extern Parameter_t HBTConfig;

inline void trim_leading_garbage(string &s, const string &garbage_list)
{
  int pos= s.find_first_not_of(garbage_list);//look for any good staff
  if( string::npos!=pos) 
	s.erase(0, pos);//s=s.substr(pos);
  else //no good staff, clear everything
	s.clear();
}
inline void trim_trailing_garbage(string &s, const string &garbage_list)
{
  int pos=s.find_first_of(garbage_list);
  if(string::npos!=pos)  
	s.erase(pos);
}
#endif
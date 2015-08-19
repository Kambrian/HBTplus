#include "config_parser.h"
void Parameter_t::SetParameterValue(const string &line)
{
  stringstream ss(line);
  string name;
  ss>>name;
//   transform(name.begin(),name.end(),name.begin(),::tolower);
#define TrySetPar(var) if(name==#var) ss>>var;
  TrySetPar(SnapshotPath)
  else TrySetPar(SnapshotFileBase)
  else TrySetPar(HaloPath)
  else TrySetPar(SubhaloPath)
  else TrySetPar(MinNumPartOfSub)
#undef TrySetPar  
  else
  {
	stringstream msg;
	msg<<"unrecognized configuration entry: "<<name<<endl;
	throw runtime_error(msg.str());
  }
}
void Parameter_t::ParseConfigFile(char* param_file)
{
  ifstream ifs;
  ifs.open(param_file);
  vector <string> lines;
  string line;
  
  cout<<"Reading configuration file "<<param_file<<endl;
  
  while(getline(ifs,line))
  {
	trim_trailing_garbage(line, "#");
	trim_leading_garbage(line, " \t");
	if(!line.empty()) SetParameterValue(line);
  }
}

Parameter_t HBTConfig;

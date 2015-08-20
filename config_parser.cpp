#include <cstdlib>
#include "config_parser.h"

void Parameter_t::SetParameterValue(const string &line)
{
  stringstream ss(line);
  string name;
  ss>>name;
//   transform(name.begin(),name.end(),name.begin(),::tolower);
#define TrySetPar(var,i) if(name==#var){ ss>>var; IsSet[i]=true; cout<<#var<<" = "<<var<<endl;}
  TrySetPar(SnapshotPath,0)
  else TrySetPar(SnapshotFileBase,1)
  else TrySetPar(HaloPath,2)
  else TrySetPar(SubhaloPath,3)
  else TrySetPar(MinNumPartOfSub,4)
  else TrySetPar(MassInMsunh,5)
  else TrySetPar(LengthInMpch,6)
  else TrySetPar(VelInKmS,7)
  else TrySetPar(BoxSize,8)
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
  CheckUnsetParameters();
  PhysicalConst::G=43.0071*(MassInMsunh/1e10)/VelInKmS/VelInKmS/LengthInMpch;
  PhysicalConst::H0=100.*(1./VelInKmS)/(1./LengthInMpch);
}
void Parameter_t::CheckUnsetParameters()
{
  for(int i=0;i<NumberOfConfigEntries;i++)
  {
	if(!IsSet[i])
	{
	  cerr<<"Error parsing configuration file: entry "<<i<<" missing\n";
	  exit(1);
	}
  }
}

Parameter_t HBTConfig;

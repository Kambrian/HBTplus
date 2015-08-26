#include <cstdlib>
#include "config_parser.h"

HBTReal PhysicalConst::G;
HBTReal PhysicalConst::H0;
Parameter_t HBTConfig;

void Parameter_t::SetParameterValue(const string &line)
{
  stringstream ss(line);
  string name;
  ss>>name;
//   transform(name.begin(),name.end(),name.begin(),::tolower);
#define TrySetPar(var,i) if(name==#var){ ss>>var; IsSet[i]=true; cout<<#var<<" = "<<var<<endl;}
  TrySetPar(SnapshotPath,0)
  else TrySetPar(HaloPath,1)
  else TrySetPar(SubhaloPath,2)
  else TrySetPar(SnapshotFileBase,3)
  else TrySetPar(MaxSnapshotIndex,4)
  else TrySetPar(BoxSize,5)
#undef TrySetPar		
#define TrySetPar(var) if(name==#var){ ss>>var; cout<<#var<<" = "<<var<<endl;}	
  else TrySetPar(MinNumPartOfSub)
  else TrySetPar(MassInMsunh)
  else TrySetPar(LengthInMpch)
  else TrySetPar(VelInKmS)
  else TrySetPar(PeriodicBoundaryOn)
  else TrySetPar(SnapshotHasIdBlock)
  else TrySetPar(ParticleIdRankStyle)
  else TrySetPar(SnapshotIdUnsigned)
#undef TrySetPar
  else if("SnapshotIdList"==name)
  {
	for(int i; ss>>i;)
	SnapshotIdList.push_back(i);
  }
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
  for(int i=0;i<NumberOfCompulsaryConfigEntries;i++)
  {
	if(!IsSet[i])
	{
	  cerr<<"Error parsing configuration file: entry "<<i<<" missing\n";
	  exit(1);
	}
  }
  if(!SnapshotIdList.empty())
  {
	int max_index=SnapshotIdList.size()-1;
	if(MaxSnapshotIndex!=max_index)
	{
	  cerr<<"Error: MaxSnapshotIndex="<<MaxSnapshotIndex<<", inconsistent with SnapshotIdList ("<<max_index+1<<" snapshots listed)\n";
	  exit(1);
	}
  }
}


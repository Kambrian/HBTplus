#ifndef CONFIG_PARSER_H_INCLUDED
#define CONFIG_PARSER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
using namespace std;

#include "datatypes.h"

namespace PhysicalConst
{//initialized after reading parameter file.
  static HBTReal G;
  static HBTReal H0;
}

#define NumberOfCompulsaryConfigEntries 6
class Parameter_t
{
public:
  string SnapshotPath;
  string HaloPath;
  string SubhaloPath;
  string SnapshotFileBase;
  int MaxSnapshotIndex;
  HBTReal BoxSize; //to check the unit of snapshot according to the BoxSize in header
  bool IsSet[NumberOfCompulsaryConfigEntries];
  
  int MinSnapshotIndex;
  int MinNumPartOfSub;
  HBTReal MassInMsunh;
  HBTReal LengthInMpch;
  HBTReal VelInKmS;
  bool PeriodicBoundaryOn;
  bool SnapshotHasIdBlock;
  bool ParticleIdRankStyle;//load particleId as id ranks
  bool SnapshotIdUnsigned;
  
  Parameter_t(): IsSet()
  {
	MinSnapshotIndex=0;
	MinNumPartOfSub=20;
	MassInMsunh=1e10;
	LengthInMpch=1e3;
	VelInKmS=1.;
	PeriodicBoundaryOn=true;
	SnapshotHasIdBlock=true;
	ParticleIdRankStyle=false;
	SnapshotIdUnsigned=false;
  }
  void ParseConfigFile(char * param_file);
  void SetParameterValue(const string &line);
  void CheckUnsetParameters();
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
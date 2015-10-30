#ifndef CONFIG_PARSER_H_INCLUDED
#define CONFIG_PARSER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include "datatypes.h"

namespace PhysicalConst
{//initialized after reading parameter file.
  extern HBTReal G;
  extern HBTReal H0;
}

#define NumberOfCompulsaryConfigEntries 7
class Parameter_t
{
public:
  /*compulsory parameters*/
  string SnapshotPath;
  string HaloPath;
  string SubhaloPath;
  string SnapshotFileBase;
  int MaxSnapshotIndex;
  HBTReal BoxSize; //to check the unit of snapshot according to the BoxSize in header
  HBTReal SofteningHalo;
  bool IsSet[NumberOfCompulsaryConfigEntries];
  
  /*optional*/
  int MinSnapshotIndex;
  int MinNumPartOfSub;
  int GroupFileVariant;
  HBTReal MassInMsunh;
  HBTReal LengthInMpch;
  HBTReal VelInKmS;
  bool PeriodicBoundaryOn;
  bool SnapshotHasIdBlock;
  bool ParticleIdRankStyle;//load particleId as id ranks
  bool ParticleIdNeedHash;//performance related; disabled if ParticleIdRankStyle is true
  bool SnapshotIdUnsigned;
  vector <int> SnapshotIdList;
  
  bool TrimNonHostParticles; //whether to trim particles outside the host halo when finding hosts
  HBTReal MajorProgenitorMassRatio; 
  HBTReal BoundMassPrecision;
  HBTReal SourceSubRelaxFactor;
  HBTReal SubCoreSizeFactor; //coresize=Nbound*CoreSizeFactor, to get center coordinates for the KineticDistance test.
  HBTInt SubCoreSizeMin; //Minimum coresize
  
  HBTReal TreeAllocFactor;
  HBTReal TreeNodeOpenAngle;
  HBTInt TreeMinNumOfCells;
  
  /*derived parameters; do not require user input*/
  HBTReal TreeNodeOpenAngleSquare;
  HBTReal TreeNodeResolution;
  HBTReal TreeNodeResolutionHalf;
  HBTReal BoxHalf; 
  
  Parameter_t(): IsSet(),SnapshotIdList()
  {
	MinSnapshotIndex=0;
	MinNumPartOfSub=20;
	GroupFileVariant=GROUP_FORMAT_GADGET3_INT;
	MassInMsunh=1e10;
	LengthInMpch=1e3;
	VelInKmS=1.;
	PeriodicBoundaryOn=true;
	SnapshotHasIdBlock=true;
	ParticleIdRankStyle=false;
	ParticleIdNeedHash=true;
	SnapshotIdUnsigned=false;
	TrimNonHostParticles=false;
	MajorProgenitorMassRatio=0.67;
	BoundMassPrecision=0.9;
	SourceSubRelaxFactor=3.;
	SubCoreSizeFactor=0.25;
	SubCoreSizeMin=20;
	TreeAllocFactor=1.; /* a value of 2 should be more than sufficient*/
	TreeNodeOpenAngle=0.45;
	TreeMinNumOfCells=500;
  }
  void ParseConfigFile(const char * param_file);
  void SetParameterValue(const string &line);
  void CheckUnsetParameters();
};

extern Parameter_t HBTConfig;
extern void ParseHBTParams(int argc, char **argv, Parameter_t &config, int &snapshot_start, int &snapshot_end);
extern void MarkHBTVersion();
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
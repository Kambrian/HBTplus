#ifndef CONFIG_PARSER_H_INCLUDED
#define CONFIG_PARSER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cmath>
#include "datatypes.h"

#define HBT_VERSION "1.15.0"

namespace PhysicalConst
{//initialized after reading parameter file.
  extern HBTReal G;
  extern HBTReal H0;
}

#define NumberOfCompulsaryConfigEntries 7
class Parameter_t
{
public:
  //remember to update SetParameterValue() and DumpParameters() accordingly if you change any parameter definition.
  /*compulsory parameters*/
  string SnapshotPath;
  string HaloPath;
  string SubhaloPath;
  string SnapshotFileBase;
  int MaxSnapshotIndex;
  HBTReal BoxSize; //to check the unit of snapshot according to the BoxSize in header
  HBTReal SofteningHalo;
  vector<bool> IsSet;
  
  /*optional*/
  string SnapshotFormat;
  string GroupFileFormat;
  int MaxConcurrentIO;
  int MinSnapshotIndex;
  int MinNumPartOfSub;
  long GroupParticleIdMask; //only used for a peculiar gadget format.
  HBTReal MassInMsunh;
  HBTReal LengthInMpch;
  HBTReal VelInKmS;
  bool PeriodicBoundaryOn;
  bool SnapshotHasIdBlock; //set to False when your snapshot is sorted according to particle id so that no id block is present.
  bool ParticleIdRankStyle;//performance related; load particleId as id ranks. not implemented yet.
  bool ParticleIdNeedHash;//performance related; disabled if ParticleIdRankStyle is true
  bool SnapshotIdUnsigned;
  bool SaveSubParticleProperties;
  bool MergeTrappedSubhalos;//whether to MergeTrappedSubhalos, see code paper for ref.
  vector <int> SnapshotIdList;
  vector <string> SnapshotNameList;
  
  HBTReal MajorProgenitorMassRatio; 
  HBTReal BoundMassPrecision;
  HBTReal SourceSubRelaxFactor;
  HBTReal SubCoreSizeFactor; //coresize=Nbound*CoreSizeFactor, to get center coordinates for the KineticDistance test.
  HBTInt SubCoreSizeMin; //Minimum coresize
  
  //only relevant if you use the FoF code shipped with HBT
  float FoFLinkParam; //dimensionless linking length b.
  int FoFSizeMin;
  float FoFDarkMatterMassFraction; 
  
  HBTReal TreeAllocFactor;
  HBTReal TreeNodeOpenAngle;
  HBTInt TreeMinNumOfCells;
  
  HBTInt MaxSampleSizeOfPotentialEstimate;
  bool RefineMostboundParticle; //whether to further improve mostbound particle accuracy in case a MaxSampleSizeOfPotentialEstimate is used. this introduces some overhead if true, but leads to more accuracy mostbound particle
  
  /*derived parameters; do not require user input*/
  HBTReal TreeNodeOpenAngleSquare;
  HBTReal TreeNodeResolution;
  HBTReal TreeNodeResolutionHalf;
  HBTReal BoxHalf; 
  bool GroupLoadedIndex; //when adding your own halo format/reader, set HBTConfig.GroupLoadedIndex=true in your read function if your group data contains particle indices rather than particle ids.
  string ConfigFile;
  
  Parameter_t(): IsSet(NumberOfCompulsaryConfigEntries,false),SnapshotIdList(),SnapshotNameList()
  {
	SnapshotFormat="gadget";
	GroupFileFormat="gadget3_int";
	MaxConcurrentIO=10;
	MinSnapshotIndex=0;
	MinNumPartOfSub=20;
	GroupParticleIdMask=0;
	MassInMsunh=1e10;
	LengthInMpch=1;
	VelInKmS=1.;
	PeriodicBoundaryOn=true;
	SnapshotHasIdBlock=true;
	ParticleIdRankStyle=false;
	ParticleIdNeedHash=true;
	SnapshotIdUnsigned=false;
	SaveSubParticleProperties=false;
#ifdef NO_STRIPPING
	MergeTrappedSubhalos=false;
#else
	MergeTrappedSubhalos=true;
#endif
	MajorProgenitorMassRatio=0.8;
	BoundMassPrecision=0.995;
	SourceSubRelaxFactor=3.;
	SubCoreSizeFactor=0.25;
	SubCoreSizeMin=20;
	FoFLinkParam=0.2;
	FoFSizeMin=20;
	FoFDarkMatterMassFraction=1.0;
	TreeAllocFactor=1.; /* a value of 2 should be more than sufficient*/
	TreeNodeOpenAngle=0.45;
	TreeMinNumOfCells=500;
	MaxSampleSizeOfPotentialEstimate=1000;//set to 0 to disable sampling
	RefineMostboundParticle=true;
	GroupLoadedIndex=false;
  }
  void ReadSnapshotNameList();
  void ParseConfigFile(const char * param_file);
  void SetParameterValue(const string &line);
  void CheckUnsetParameters();
  void DumpParameters(const string &path);
  void DumpParameters();
};

extern Parameter_t HBTConfig;
extern void ParseHBTParams(int argc, char **argv, Parameter_t &config, int &snapshot_start, int &snapshot_end);
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

#define NEAREST(x) (((x)>HBTConfig.BoxHalf)?((x)-HBTConfig.BoxSize):(((x)<-HBTConfig.BoxHalf)?((x)+HBTConfig.BoxSize):(x)))
inline HBTReal PeriodicDistance(const HBTxyz &x, const HBTxyz &y)
{
	HBTxyz dx;
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];
	dx[2]=x[2]-y[2];
	if(HBTConfig.PeriodicBoundaryOn)
	{
	  dx[0]=NEAREST(dx[0]);
	  dx[1]=NEAREST(dx[1]);
	  dx[2]=NEAREST(dx[2]);
	}
	return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
}
#endif
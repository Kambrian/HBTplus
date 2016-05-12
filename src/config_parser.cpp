#include <cstdlib>
#include "config_parser.h"
#include "mymath.h"

namespace PhysicalConst
{
HBTReal G;
HBTReal H0;
}

Parameter_t HBTConfig;

void Parameter_t::SetParameterValue(const string &line)
{
  stringstream ss(line);
  string name;
  ss>>name;
//   transform(name.begin(),name.end(),name.begin(),::tolower);
#define TrySetPar(var,i) if(name==#var){ ss>>var; IsSet[i]=true;}
  TrySetPar(SnapshotPath,0)
  else TrySetPar(HaloPath,1)
  else TrySetPar(SubhaloPath,2)
  else TrySetPar(SnapshotFileBase,3)
  else TrySetPar(MaxSnapshotIndex,4)
  else TrySetPar(BoxSize,5)
  else TrySetPar(SofteningHalo,6)
#undef TrySetPar		
#define TrySetPar(var) if(name==#var) ss>>var;
  else TrySetPar(SnapshotFormat)
  else TrySetPar(MaxConcurrentIO)
  else TrySetPar(MinSnapshotIndex)
  else TrySetPar(MinNumPartOfSub)
  else TrySetPar(GroupFileVariant)
//   else TrySetPar(GroupParticleIdMask)
  else if(name=="GroupParticleIdMask")
  {
        ss>>hex>>GroupParticleIdMask>>dec;
        cout<<"GroupParticleIdMask = "<<hex<<GroupParticleIdMask<<dec<<endl;
  }	
  else TrySetPar(MassInMsunh)
  else TrySetPar(LengthInMpch)
  else TrySetPar(VelInKmS)
  else TrySetPar(PeriodicBoundaryOn)
  else TrySetPar(SnapshotHasIdBlock)
  else TrySetPar(ParticleIdRankStyle)
  else TrySetPar(ParticleIdNeedHash)
  else TrySetPar(SnapshotIdUnsigned)
  else TrySetPar(SaveSubParticleProperties)
  else TrySetPar(MajorProgenitorMassRatio)
#ifdef ALLOW_BINARY_SYSTEM
  else TrySetPar(BinaryMassRatioLimit)
#endif
  else TrySetPar(BoundMassPrecision)
  else TrySetPar(SourceSubRelaxFactor)
  else TrySetPar(SubCoreSizeFactor)
  else TrySetPar(SubCoreSizeMin)
  else TrySetPar(TreeAllocFactor)
  else TrySetPar(TreeNodeOpenAngle)
  else TrySetPar(TreeMinNumOfCells)
  else TrySetPar(MaxSampleSizeOfPotentialEstimate)
  else TrySetPar(RefineMostboundParticle)
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

void Parameter_t::ParseConfigFile(const char * param_file)
{
  ifstream ifs;
  ifs.open(param_file);
  if(!ifs.is_open())
  {
	cerr<<"Error opening config file "<<param_file<<endl;
	exit(1);
  }
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
  
  if(ParticleIdRankStyle) ParticleIdNeedHash=false;
  
  BoxHalf=BoxSize/2.;
  TreeNodeResolution=SofteningHalo*0.1;
  TreeNodeResolutionHalf=TreeNodeResolution/2.;
  TreeNodeOpenAngleSquare=TreeNodeOpenAngle*TreeNodeOpenAngle;
  
  ReadSnapshotNameList();
}
void Parameter_t::ReadSnapshotNameList()
{//to specify snapshotnamelist, create a file "snapshotlist.txt" under SubhaloPath, and list the filenames inside, one per line.
  string snaplist_filename=SubhaloPath+"/snapshotlist.txt";
  ifstream ifs;
  ifs.open(snaplist_filename);
  if(ifs.is_open())
  {
	cout<<"Found SnapshotNameList file "<<snaplist_filename<<endl;
	
	string line;	
	while(getline(ifs,line))
	{
	  trim_trailing_garbage(line, "#");
	  istringstream ss(line);
	  string name;
	  ss>>name;
	  if(!name.empty())
	  {
		cout<<name<<endl;
		SnapshotNameList.push_back(name);
	  }
	}
  }
  if(SnapshotNameList.size())
	assert(SnapshotNameList.size()==MaxSnapshotIndex+1);
}

void Parameter_t::CheckUnsetParameters()
{
  for(int i=0;i<IsSet.size();i++)
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

void ParseHBTParams(int argc, char **argv, Parameter_t &config, int &snapshot_start, int &snapshot_end)
{
  if(argc<2)
  {
	cerr<<"Usage: "<<argv[0]<<" [param_file] <snapshot_start> <snapshot_end>\n";
	exit(1);
  }
  config.ParseConfigFile(argv[1]);
  if(2==argc)
  {
	snapshot_start=config.MinSnapshotIndex;
	snapshot_end=config.MaxSnapshotIndex;
  }
  else
  {
  snapshot_start=atoi(argv[2]);
  if(argc>3)
	snapshot_end=atoi(argv[3]);
  else
	snapshot_end=snapshot_start;
  }
  cout<<"Running "<<argv[0]<<" from snapshot "<<snapshot_start<<" to "<<snapshot_end<<" using configuration file "<<argv[1]<<endl;
}

void Parameter_t::DumpParameters()
{
#ifndef HBT_VERSION
  #define HBT_VERSION "unknown"
  cout<<"Warning: HBT_VERSION not set. Better write down which version you are using.\n";
#endif
  string filename=SubhaloPath+"/VER"+HBT_VERSION+".param";
  ofstream version_file(filename, ios::out|ios::trunc);
  if(!version_file.is_open())
  {
	cerr<<"Error opening "<<filename<<" for parameter dump."<<endl;
	exit(1);
  }
#define DumpPar(var) version_file<<#var<<"  "<<var<<endl;
  DumpPar(SnapshotPath)
  DumpPar(HaloPath)
  DumpPar(SubhaloPath)
  DumpPar(SnapshotFileBase)
  DumpPar(MaxSnapshotIndex)
  DumpPar(BoxSize)
  DumpPar(SofteningHalo)
  
  /*optional*/
  DumpPar(SnapshotFormat)
  DumpPar(MaxConcurrentIO)
  DumpPar(MinSnapshotIndex)
  DumpPar(MinNumPartOfSub)
  DumpPar(GroupFileVariant)
  if(GroupParticleIdMask)
	version_file<<"GroupParticleIdMask "<<hex<<GroupParticleIdMask<<dec<<endl;
  DumpPar(MassInMsunh)
  DumpPar(LengthInMpch)
  DumpPar(VelInKmS)
  DumpPar(PeriodicBoundaryOn)
  DumpPar(SnapshotHasIdBlock)
  DumpPar(ParticleIdRankStyle)
  DumpPar(ParticleIdNeedHash)
  DumpPar(SnapshotIdUnsigned)
  DumpPar(SaveSubParticleProperties)
  if(SnapshotIdList.size())
  {
  version_file<<"SnapshotIdList";
  for(auto && i: SnapshotIdList)
	version_file<<" "<<i;
  version_file<<endl;
  }
  if(SnapshotNameList.size())
  {
  version_file<<"#SnapshotNameList";
  for(auto && i: SnapshotNameList)
	version_file<<" "<<i;
  version_file<<endl;
  }
  
  DumpPar(MajorProgenitorMassRatio) 
#ifdef ALLOW_BINARY_SYSTEM
  DumpPar(BinaryMassRatioLimit)
#endif
  DumpPar(BoundMassPrecision)
  DumpPar(SourceSubRelaxFactor)
  DumpPar(SubCoreSizeFactor) 
  DumpPar(SubCoreSizeMin)
  
  DumpPar(TreeAllocFactor)
  DumpPar(TreeNodeOpenAngle)
  DumpPar(TreeMinNumOfCells)
  
  DumpPar(MaxSampleSizeOfPotentialEstimate)
  DumpPar(RefineMostboundParticle)
#undef DumpPar  
  version_file.close();
}

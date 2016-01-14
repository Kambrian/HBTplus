#include <cstdlib>
#include "config_parser.h"

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
//   else TrySetPar(SnapshotNoMassBlock)
  else TrySetPar(ParticleIdRankStyle)
  else TrySetPar(ParticleIdNeedHash)
  else TrySetPar(SnapshotIdUnsigned)
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
  if(!ifs.is_open())//or ifs.fail()
  {
	cerr<<"Error: failed to open configuration: "<<param_file<<endl;
	exit(1);
  }
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
  
  if(ParticleIdRankStyle) ParticleIdNeedHash=false;
  
  BoxHalf=BoxSize/2.;
  TreeNodeResolution=SofteningHalo*0.1;
  TreeNodeResolutionHalf=TreeNodeResolution/2.;
  TreeNodeOpenAngleSquare=TreeNodeOpenAngle*TreeNodeOpenAngle;
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

void Parameter_t::BroadCast(MpiWorker_t &world, int root)
/*sync parameters and physical consts across*/
{
  #define _SyncVec(x,t) world.SyncContainer(x,t,root)	
  #define _SyncAtom(x,t) world.SyncAtom(x,t,root)
  #define _SyncBool(x) world.SyncAtomBool(x, root)
  #define _SyncVecBool(x) world.SyncVectorBool(x, root)
  #define _SyncReal(x) _SyncAtom(x, MPI_HBT_REAL)
  
  _SyncVec(SnapshotPath, MPI_CHAR);
  _SyncVec(HaloPath, MPI_CHAR);
  _SyncVec(SubhaloPath, MPI_CHAR);
  _SyncVec(SnapshotFileBase, MPI_CHAR);
  _SyncAtom(MaxSnapshotIndex, MPI_INT);
  _SyncReal(BoxSize);
  _SyncReal(SofteningHalo);
  _SyncVecBool(IsSet);
  
  _SyncAtom(MinSnapshotIndex, MPI_INT);
  _SyncAtom(MinNumPartOfSub, MPI_INT);
  _SyncAtom(GroupFileVariant, MPI_INT);
  _SyncAtom(GroupParticleIdMask, MPI_LONG);
  _SyncReal(MassInMsunh);
  _SyncReal(LengthInMpch);
  _SyncReal(VelInKmS);
  _SyncBool(PeriodicBoundaryOn);
  _SyncBool(SnapshotHasIdBlock);
  _SyncBool(ParticleIdRankStyle);
  _SyncBool(ParticleIdNeedHash);
  _SyncBool(SnapshotIdUnsigned);
  _SyncVec(SnapshotIdList, MPI_INT);
  
  _SyncReal(MajorProgenitorMassRatio);
  #ifdef ALLOW_BINARY_SYSTEM
  _SyncReal(BinaryMassRatioLimit);
  #endif
  _SyncReal(BoundMassPrecision);
  _SyncReal(SourceSubRelaxFactor);
  _SyncReal(SubCoreSizeFactor);
  _SyncAtom(SubCoreSizeMin, MPI_HBT_INT);
  
  _SyncReal(TreeAllocFactor);
  _SyncReal(TreeNodeOpenAngle);
  _SyncAtom(TreeMinNumOfCells, MPI_HBT_INT);
  
  _SyncReal(TreeNodeOpenAngleSquare);
  _SyncReal(TreeNodeResolution);
  _SyncReal(TreeNodeResolutionHalf);
  _SyncReal(BoxHalf);
  
  //---------------end sync params-------------------------//	
  
  _SyncReal(PhysicalConst::G);
  _SyncReal(PhysicalConst::H0);
  
  #undef _SyncVec
  #undef _SyncAtom 
  #undef _SyncBool 
  #undef _SyncVecBool 
  #undef _SyncReal 
}
void Parameter_t::DumpParameters()
{
#ifndef HBT_VERSION
	cout<<"Warning: HBT_VERSION unknown.\n";
#define HBT_VERSION "unknown"
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
  if(SnapshotIdList.size())
  {
  version_file<<"SnapshotIdList";
  for(auto && i: SnapshotIdList)
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
#undef DumpPar  
  version_file.close();
}

using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <chrono>

#include "../snapshot.h"
#include "../mymath.h"
#include "../boost_mpi.h"


#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
#define SkipPositionBlock SkipBlock
#define SkipVelocityBlock SkipBlock	 
#define SkipIdBlock SkipBlock
#define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)

#define BLOCK_ID_POS 0
#define BLOCK_ID_VEL 1
#define BLOCK_ID_PID 2
#define BLOCK_ID_MASS 3
static int BLOCK_DIMENSION[]={3,3,1,1};

inline size_t ParticleSnapshot_t::SkipBlock(FILE *fp)
{
  int blocksize,blocksize2;
  
  ReadBlockSize(blocksize);	  
  fseek(fp, blocksize, SEEK_CUR);
  ReadBlockSize(blocksize2);
  assert(blocksize==blocksize2);
  return blocksize;
}
void ParticleSnapshot_t::GetFileName(int ifile, string &filename)
{
  FILE *fp;
  char buf[1024];
  
  sprintf(buf,"%s/snapdir_%03d/%s_%03d.%d",HBTConfig.SnapshotPath.c_str(),SnapshotId,HBTConfig.SnapshotFileBase.c_str(),SnapshotId,ifile);
  if(ifile==0)
	if(!file_exist(buf)) sprintf(buf,"%s/%s_%03d",HBTConfig.SnapshotPath.c_str(),HBTConfig.SnapshotFileBase.c_str(),SnapshotId); //try the other convention
  if(!file_exist(buf)) sprintf(buf,"%s/%s_%03d.%d",HBTConfig.SnapshotPath.c_str(),HBTConfig.SnapshotFileBase.c_str(),SnapshotId,ifile); //try the other convention
  if(!file_exist(buf)) sprintf(buf,"%s/%d/%s.%d",HBTConfig.SnapshotPath.c_str(),SnapshotId,HBTConfig.SnapshotFileBase.c_str(),ifile);//for BJL's RAMSES output
  if(!file_exist(buf))
  {
	cerr<<"Failed to find a snapshot file at "<<buf<<endl;
	exit(1);
  }
  filename=buf;
}
bool ParticleSnapshot_t::ReadFileHeader(FILE *fp, SnapshotHeader_t &header)
{
  //read the header part, assign header extensions, and check byteorder
  int headersize(SNAPSHOT_HEADER_SIZE),headersize_byteswap(SNAPSHOT_HEADER_SIZE);
  swap_Nbyte(&headersize_byteswap,1,sizeof(headersize_byteswap));

  bool NeedByteSwap;
  int dummy,dummy2;
  size_t tmp_size=fread(&dummy,sizeof(dummy),1,fp);
  if(dummy==headersize)
	NeedByteSwap=false;
  else if(dummy==headersize_byteswap)
	NeedByteSwap=true;
  else
  {
	cerr<<"endianness check failed for header\n file format not expected:"<<dummy<<" not match headersize "<<headersize<<" or "<<headersize_byteswap<<endl<<flush;
	exit(1);
  }
  dummy=headersize;
  
  myfread(header.npart,sizeof(int),NUMBER_OF_PARTICLE_TYPES,fp);
  myfread(header.mass,sizeof(double),NUMBER_OF_PARTICLE_TYPES,fp);
  myfread(&header.ScaleFactor,sizeof(double),1,fp);
  myfread(&header.redshift,sizeof(double),1,fp);
  myfread(&header.flag_sfr,sizeof(int),1,fp);
  myfread(&header.flag_feedback,sizeof(int),1,fp);
  myfread(header.npartTotal,sizeof(int),NUMBER_OF_PARTICLE_TYPES,fp);
  myfread(&header.flag_cooling,sizeof(int),1,fp);
  myfread(&header.num_files,sizeof(int),1,fp);
  myfread(&header.BoxSize,sizeof(double),1,fp);
  myfread(&header.Omega0,sizeof(double),1,fp);
  myfread(&header.OmegaLambda,sizeof(double),1,fp);
  myfread(&header.HubbleParam,sizeof(double),1,fp);
  fseek(fp,headersize+sizeof(int),SEEK_SET);
  myfread(&dummy2,sizeof(dummy2),1,fp);
  if(dummy!=dummy2)
  {
	cerr<<"Error: record brackets not match for header!"<<dummy<<dummy2<<endl;
	exit(1);
  } 
  
  return NeedByteSwap;
}
HBTInt ParticleSnapshot_t::ReadNumberOfDMParticles(int ifile)
{
  FILE * fp;
  string filename;
  GetFileName( ifile, filename);
  myfopen(fp, filename.c_str(), "r");
  SnapshotHeader_t header;
  ReadFileHeader(fp, header);
  fclose(fp);
  return header.npart[1];
}
void ParticleSnapshot_t::LoadHeader(int ifile)
{
  //read the header part, assign header extensions, and check byteorder

  FILE *fp; 
  string filename;
  GetFileName( ifile, filename);
  
  myfopen(fp, filename.c_str(), "r");
  NeedByteSwap=ReadFileHeader(fp, Header);
  
  if(Header.BoxSize!=HBTConfig.BoxSize)
  {
	cerr<<"BoxSize not match input: read "<<Header.BoxSize<<"; expect "<<HBTConfig.BoxSize<<endl;
	cerr<<"Maybe the length unit differ? Excpected unit: "<<HBTConfig.LengthInMpch<<" Msol/h\n";
	exit(1);
  }
  
  //set datatypes
  HBTInt NumPartInFile=0;
  for(int i=0;i<NUMBER_OF_PARTICLE_TYPES;i++) 
	NumPartInFile+=Header.npart[i];
  int blocksize=SkipPositionBlock(fp);
  RealTypeSize=blocksize/NumPartInFile/3;
  assert(sizeof(float)==RealTypeSize||sizeof(double)==RealTypeSize);
  blocksize=SkipVelocityBlock(fp);
  assert(blocksize==RealTypeSize*NumPartInFile*3);
  if(sizeof(HBTReal)<RealTypeSize)
	cerr<<"WARNING: loading size "<<RealTypeSize<<" float in snapshot with size "<<sizeof(HBTReal)<<" float in HBT. possible loss of accuracy.\n Please use ./HBTdouble unless you know what you are doing.";
  
  if(HBTConfig.SnapshotHasIdBlock)
  {
	blocksize=SkipIdBlock(fp);
	IntTypeSize=blocksize/NumPartInFile;
	assert(sizeof(int)==IntTypeSize||sizeof(long)==IntTypeSize);
	assert(sizeof(HBTInt)>=IntTypeSize);
	if(sizeof(HBTInt)<IntTypeSize)
	  cerr<<"WARNING: loading size "<<IntTypeSize<<" integer in snapshot with size "<<sizeof(HBTInt)<<" int in HBT. possible data overflow.\n Please use ./HBTdouble unless you know what you are doing.";
  }
  
  fclose(fp);
  
  SetEpoch(Header.ScaleFactor, Header.Omega0, Header.OmegaLambda);
  
  //npartTotal is not reliable
  NumberOfParticles=0;
  for(int iFile=0;iFile<Header.num_files;iFile++)
  {
	int n=ReadNumberOfDMParticles( iFile);
	NumberOfDMParticleInFiles.push_back(n);
	OffsetOfDMParticleInFiles.push_back(NumberOfParticles);
	NumberOfParticles+=n;
  }
  
}

void ParticleSnapshot_t::Load(int snapshot_index, bool fill_particle_hash)
{ 
  SetSnapshotIndex(snapshot_index);
  LoadHeader();
  
  mpi::communicator world;
  
  int nfiles_read=round(1.*Header.num_files/world.size());
  int nfiles_skip=nfiles_read*world.rank();
  int nfiles_end=min(nfiles_read+nfiles_skip, Header.num_files);
 
  accumulate(NumberOfDMParticleInFiles.begin()+nfiles_read, NumberOfDMParticleInFiles.begin()+nfiles_end, NumberOfParticles);
  if(LoadFlag.Id)  
	ParticleId.reserve(NumberOfParticles);
  if(LoadFlag.Pos)
	ComovingPosition.reserve(NumberOfParticles);
  if(LoadFlag.Vel)
	PhysicalVelocity.reserve(NumberOfParticles);
  if(LoadFlag.Mass&&0.==Header.mass[1]) 
	ParticleMass.reserve(NumberOfParticles);
  
  for(int iFile=nfiles_skip; (iFile<nfiles_read+nfiles_skip)&&(iFile<Header.num_files); iFile++)
  {
	ReadFile(iFile);
  }
 
 if(LoadFlag.Id&&fill_particle_hash)
   FillParticleHash();
   
 if(HBTConfig.PeriodicBoundaryOn)//regularize coord
  {
	for(HBTInt i=0;i<NumberOfParticles;i++)
	  for(int j=0;j<3;j++)
		ComovingPosition[i][j]=position_modulus(ComovingPosition[i][j], HBTConfig.BoxSize);
  }
  HBTReal velocity_scale=sqrt(Header.ScaleFactor);
  for(HBTInt i=0;i<NumberOfParticles;i++)
	for(int j=0;j<3;j++)
	  PhysicalVelocity[i][j]*=velocity_scale;

/*
  if((NumberOfParticles>1)&&(ParticleMass[0]!=ParticleMass[1]))
  {
	cout<<"Error: DM particles have different mass? not supported!\n";
	exit(1);
  }
  else
	Header.mass[1]=ParticleMass[0];*/

  cout<<"Finished Reading on thread "<<world.rank()<<endl;
  cout<<GetNumberOfParticles()<<" Particles loaded\n";
  cout<<GetParticleId(10)<<endl;
  cout<<GetComovingPosition(10)<<endl;
  cout<<GetParticleMass(10)<<','<<GetParticleMass(100)<<endl;
  cout<<GetParticleIndex(GetParticleId(10))<<endl;
}
template <class T, class U>
void ParticleSnapshot_t::ReadScalarBlock(FILE *fp, size_t n_read, size_t n_skip, vector <U> &x)
{
	FortranBlock <T> block(fp, n_read, n_skip, NeedByteSwap);
	x.assign(block.begin(), block.end());
}
template <class T>
void ParticleSnapshot_t::ReadXyzBlock(FILE *fp, size_t n_read, size_t n_skip, vector <HBTxyz> &x)
{
	FortranBlock <T> block(fp, n_read*3, n_skip*3, NeedByteSwap);
	//lacking assign operator for HBTxyz, has to manually assign
	HBTInt n_old=x.size();
	x.resize(n_old+n_read);
	auto p=block.data_reshape();
	for(auto it=x.begin()+n_read;it<x.end();it++)
	{
	  for(int j=0;j<3;j++)
		(*it)[j]=(*p)[j];
	  p++;
	}
}

void ParticleSnapshot_t::ReadFile(int iFile)
{ 
  FILE *fp; 
  string filename;
  GetFileName( iFile, filename);
  myfopen(fp, filename.c_str(), "r");
  SnapshotHeader_t header;
  ReadFileHeader(fp, header);
  size_t n_read=header.npart[1], n_skip=header.npart[0];
  
  if(LoadFlag.Pos)
  {
	if(RealTypeSize==4)
	  ReadXyzBlock<float>(fp, n_read, n_skip, ComovingPosition);
	else
	  ReadXyzBlock<double>(fp, n_read, n_skip, ComovingPosition);
  }
  else
	SkipBlock(fp);
  
  if(LoadFlag.Vel)
  {
	if(RealTypeSize==4)
	  ReadXyzBlock<float>(fp, n_read, n_skip, PhysicalVelocity);
	else
	  ReadXyzBlock<double>(fp, n_read, n_skip, PhysicalVelocity);
  }
  else
	SkipBlock(fp);
  
  if(LoadFlag.Id&&HBTConfig.SnapshotHasIdBlock)
  {
	if(HBTConfig.ParticleIdRankStyle)
	{
	  cout<<"Error: ParticleIdRankStyle not implemented yet\n";
	  exit(1);
	}
	
	if(IntTypeSize==4)
	{
	  if(HBTConfig.SnapshotIdUnsigned)//unsigned int
		ReadScalarBlock<unsigned>(fp, n_read, n_skip, ParticleId);
	  else
		ReadScalarBlock<int>(fp, n_read, n_skip, ParticleId);
	}
	else
	  ReadScalarBlock<long>(fp, n_read, n_skip, ParticleId);
  }
  else
	SkipBlock(fp);
  
  #define MassDataPresent(i) (0==header.mass[i])
  if(LoadFlag.Mass&&MassDataPresent(1))
  {
	if(!MassDataPresent(0)) n_skip=0;
	
	if(RealTypeSize==4)
	  ReadScalarBlock<float>(fp, n_read, n_skip, ParticleMass);
	else
	  ReadScalarBlock<double>(fp, n_read, n_skip, ParticleMass);
  }
  else
  {
	for(int i=0;i<NUMBER_OF_PARTICLE_TYPES;i++)
	  if(MassDataPresent(i))
	  {
		SkipBlock(fp);
		break;
	  }
  }
  #undef MassDataPresent
  
  if(feof(fp))
  {
	cout<<"Error: End-of-File when reading "<<filename<<endl;
	exit(1);
  }
  
  fclose(fp);
}

void ParticleSnapshot_t::SetLoadFlags(bool load_id, bool load_pos, bool load_vel, bool load_mass)
/* set the flags to only load some blocks; only the blocks with a true flag will be loaded. */
{
  LoadFlag.Id=load_id;
  LoadFlag.Pos=load_pos;
  LoadFlag.Vel=load_vel;
  LoadFlag.Mass=load_mass;
}

void ParticleSnapshot_t::Clear()
/*reset to empty*/
{
#define RESET(x, T) {vector <T>().swap(x);}
  RESET(ParticleId, HBTInt);
  RESET(ComovingPosition, HBTxyz);
  RESET(PhysicalVelocity, HBTxyz);
  RESET(ParticleMass, HBTReal);
#undef RESET 
  ClearParticleHash();//even if you don't do this, the destructor will still clean up the memory.
//   cout<<NumberOfParticles<<" particles cleared from snapshot "<<SnapshotIndex<<endl;
  NumberOfParticles=0;
}

#ifdef TEST_snapshot_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
    mpi::environment env;
  mpi::communicator world;
#ifdef _OPENMP
 omp_set_nested(0);
#endif
   
  int snapshot_start, snapshot_end;
  if(0==world.rank())
  {
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
	mkdir(HBTConfig.SubhaloPath.c_str(), 0755);
	MarkHBTVersion();
  }
  broadcast(world, HBTConfig, 0);
  
  ParticleSnapshot_t snapshot;
  snapshot.Load(HBTConfig.MaxSnapshotIndex, true);

  cout<<"Finished Loading on thread "<<world.rank()<<endl;
  cout<<snapshot.GetNumberOfParticles()<<endl;
  cout<<snapshot.GetParticleId(10)<<endl;
  cout<<snapshot.GetComovingPosition(10)<<endl;
  cout<<snapshot.GetParticleMass(10)<<','<<snapshot.GetParticleMass(100)<<endl;
  cout<<snapshot.GetParticleIndex(snapshot.GetParticleId(10))<<endl<<endl;
  
  return 0;
}
#endif
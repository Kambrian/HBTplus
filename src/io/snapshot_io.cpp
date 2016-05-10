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


#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
#define SkipPositionBlock SkipBlock
#define SkipVelocityBlock SkipBlock	 
#define SkipIdBlock SkipBlock
#define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)
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
  myfread(&header.OmegaM0,sizeof(double),1,fp);
  myfread(&header.OmegaLambda0,sizeof(double),1,fp);
  myfread(&header.HubbleParam,sizeof(double),1,fp);
  fseek(fp,headersize+sizeof(int),SEEK_SET);
  myfread(&dummy2,sizeof(dummy2),1,fp);
  if(dummy!=dummy2)
  {
	fprintf(stderr,"error!record brackets not match for header!\t%d,%d\n",dummy,dummy2);
	exit(1);
  } 
  
#ifdef MAJOR_MERGER_PATCH
/*work around for the buggy non-conforming major-merger snapshot*/ 
  header.num_files=1;
  header.BoxSize=250.;
#endif  
  
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
  
  SetEpoch(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);
  
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
  Clear();
  SetSnapshotIndex(snapshot_index);
  
  LoadHeader();
  if(LoadFlag.Id)
  {
	LoadId();
	if(fill_particle_hash)
	  FillParticleHash();
  }
  if(LoadFlag.Pos)
	LoadPosition();
  if(LoadFlag.Vel)
	LoadVelocity();
  if(LoadFlag.Mass)
	if(0.==Header.mass[1]) LoadMass();
}

void ParticleSnapshot_t::SetLoadFlags(bool load_id, bool load_pos, bool load_vel, bool load_mass)
/* set the flags to only load some blocks; only the blocks with a true flag will be loaded. */
{
  LoadFlag.Id=load_id;
  LoadFlag.Pos=load_pos;
  LoadFlag.Vel=load_vel;
  LoadFlag.Mass=load_mass;
}

size_t ParticleSnapshot_t::ReadBlock(FILE *fp, void *block, const size_t n_read, const size_t n_skip_before, const size_t n_skip_after)
{//read n_read members from the current block of fp into block. skip n_skip_* before and after. 
  //return member size in block. 
 //has to specify the number of members in order to know the member size for byteswap
	  int blocksize,blocksize2;
	  
	  ReadBlockSize(blocksize);	  
	  size_t block_member_size=blocksize/(n_read+n_skip_before+n_skip_after);
	  assert(4==block_member_size||8==block_member_size);
	  fseek(fp, n_skip_before*block_member_size, SEEK_CUR);
	  myfread(block, block_member_size, n_read, fp);
	  fseek(fp, n_skip_after*block_member_size, SEEK_CUR);
	  
	  ReadBlockSize(blocksize2);
	  assert(blocksize==blocksize2);
	  
	  return block_member_size;
}
void * ParticleSnapshot_t::LoadBlock( int block_id, size_t element_size, int dimension, bool is_massblock)
{
  char * buf=static_cast<char *>(::operator new(element_size*NumberOfParticles*dimension));
  //#pragma omp parallel //can do task parallelization here.
  for(int iFile=0; iFile<Header.num_files; iFile++)
  {
	FILE *fp; 
	string filename;
	GetFileName( iFile, filename);
	myfopen(fp, filename.c_str(), "r");
	
	SnapshotHeader_t header;
	ReadFileHeader(fp, header);
	for(int iblock=0;iblock<block_id;iblock++)
	  SkipBlock(fp);
	size_t n_read=header.npart[1], n_skip_before=0, n_skip_after=0;
#define MassDataPresent(i) (0==header.mass[i])
	if(!is_massblock||MassDataPresent(0))
	  n_skip_before=header.npart[0];
	for(int i=2;i<NUMBER_OF_PARTICLE_TYPES;i++)
	{
	 if(!is_massblock||MassDataPresent(i))
		n_skip_after+=header.npart[i];
	}
#undef MassDataPresent
	size_t buf_member_size=ReadBlock(fp, buf+OffsetOfDMParticleInFiles[iFile]*dimension*element_size, n_read*dimension, n_skip_before*dimension, n_skip_after*dimension);
	assert(buf_member_size==element_size);
	
	if(feof(fp))
	{
	  cerr<<"error:End-of-File in "<<filename<<endl;
	  exit(1);  
	}
	
	fclose(fp);
  }
  
  return static_cast<void *>(buf);
}

template <class T, class U>
static void AssignBlock(T *dest, const U *source, const size_t n)
{
  for(size_t i=0;i<n;i++)
  {
	dest[i]=source[i];
  }
}

#define BLOCK_ID_POS 0
#define BLOCK_ID_VEL 1
#define BLOCK_ID_PID 2
#define BLOCK_ID_MASS 3

void ParticleSnapshot_t::LoadId()
{//only do this after LoadHeader().
  
  if(!HBTConfig.SnapshotHasIdBlock) return;
  
  void * buf=LoadBlock( BLOCK_ID_PID, IntTypeSize, 1);
  
  if(HBTConfig.ParticleIdRankStyle)
  {
	cerr<<"Error: ParticleIdRankStyle not implemented yet\n";
	exit(1);
  }
  else
  {
	if(HBTConfig.SnapshotIdUnsigned)//unsigned int
	  AssignBlock(ParticleId, static_cast<unsigned *>(buf), NumberOfParticles);
	else
	{
	  if(sizeof(HBTInt)==IntTypeSize)
		ParticleId=static_cast<HBTInt *>(buf);
	  else
	  {
		ParticleId=static_cast<HBTInt *>(::operator new(sizeof(HBTInt)*NumberOfParticles));
		if(sizeof(int)==IntTypeSize)
		  AssignBlock(ParticleId, static_cast<int *>(buf), NumberOfParticles);
		else
		  AssignBlock(ParticleId, static_cast<long *>(buf), NumberOfParticles);
		::operator delete(buf);
	  }
	}
  }
}

void ParticleSnapshot_t::LoadPosition()
{
  void * buf=LoadBlock( BLOCK_ID_POS, RealTypeSize, 3);
  
  if(sizeof(HBTReal)==RealTypeSize)
	ComovingPosition=static_cast<HBTxyz *>(buf);
  else
  {
	ComovingPosition=static_cast<HBTxyz *>(::operator new(sizeof(HBTxyz)*NumberOfParticles));
	if(sizeof(float)==RealTypeSize)
	  AssignBlock(reinterpret_cast<HBTReal *>(ComovingPosition), static_cast<float *>(buf), NumberOfParticles*3);
	else
	  AssignBlock(reinterpret_cast<HBTReal *>(ComovingPosition), static_cast<double *>(buf), NumberOfParticles*3);
	::operator delete(buf);
  }
  
  if(HBTConfig.PeriodicBoundaryOn)//regularize coord
  {
	for(HBTInt i=0;i<NumberOfParticles;i++)
	  for(int j=0;j<3;j++)
		ComovingPosition[i][j]=position_modulus(ComovingPosition[i][j], HBTConfig.BoxSize);
  }
}

void ParticleSnapshot_t::LoadVelocity()
{
  void * buf=LoadBlock( BLOCK_ID_VEL, RealTypeSize, 3);
  
  if(sizeof(HBTReal)==RealTypeSize)
	PhysicalVelocity=static_cast<HBTxyz *>(buf);
  else
  {
	PhysicalVelocity=static_cast<HBTxyz *>(::operator new(sizeof(HBTxyz)*NumberOfParticles));
	if(sizeof(float)==RealTypeSize)
	  AssignBlock(reinterpret_cast<HBTReal *>(PhysicalVelocity), static_cast<float *>(buf), NumberOfParticles*3);
	else
	  AssignBlock(reinterpret_cast<HBTReal *>(PhysicalVelocity), static_cast<double *>(buf), NumberOfParticles*3);
	::operator delete(buf);
  }
  
  HBTReal velocity_scale=sqrt(Header.ScaleFactor);
  for(HBTInt i=0;i<NumberOfParticles;i++)
	for(int j=0;j<3;j++)
	  PhysicalVelocity[i][j]*=velocity_scale;
}

void ParticleSnapshot_t::LoadMass()
{
  void * buf=LoadBlock( BLOCK_ID_MASS, RealTypeSize, 1, true);
  
  if(sizeof(HBTReal)==RealTypeSize)
	ParticleMass=static_cast<HBTReal *>(buf);
  else
  {
	ParticleMass=static_cast<HBTReal *>(::operator new(sizeof(HBTReal)*NumberOfParticles));
	if(sizeof(float)==RealTypeSize)
	  AssignBlock(ParticleMass, static_cast<float *>(buf), NumberOfParticles);
	else
	  AssignBlock(ParticleMass, static_cast<double *>(buf), NumberOfParticles);
	::operator delete(buf);
  }/*
  if((NumberOfParticles>1)&&(ParticleMass[0]!=ParticleMass[1]))
  {
	cout<<"Error: DM particles have different mass? not supported!\n";
	exit(1);
  }
  else
	Header.mass[1]=ParticleMass[0];*/
}

void ParticleSnapshot_t::Clear()
/*reset to empty*/
{
#define RESET(x) {::operator delete(x);x=nullptr;}
  RESET(ParticleId);
  RESET(ComovingPosition);
  RESET(PhysicalVelocity);
  RESET(ParticleMass);
#undef RESET 
  ClearParticleHash();//even if you don't do this, the destructor will still clean up the memory.
//   cout<<NumberOfParticles<<" particles cleared from snapshot "<<SnapshotIndex<<endl;
  NumberOfParticles=0;
}

#ifdef TEST_snapshot_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  HBTConfig.ParseConfigFile(argv[1]);
  ParticleSnapshot_t snapshot;
  snapshot.Load(HBTConfig.MaxSnapshotIndex, true, true, false, true, true);
  cout<<snapshot.GetNumberOfParticles()<<endl;
  cout<<snapshot.GetParticleId(10)<<endl;
  cout<<snapshot.GetComovingPosition(10)<<endl;
  cout<<snapshot.GetParticleMass(10)<<','<<snapshot.GetParticleMass(100)<<endl;
  cout<<snapshot.GetParticleIndex(snapshot.GetParticleId(10))<<endl;
  return 0;
}
#endif
using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "snapshot.h"
#include "../mymath.h"


#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
#define SkipPositionBlock SkipBlock
#define SkipVelocityBlock SkipBlock	 
#define SkipIdBlock SkipBlock
#define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)
inline size_t Snapshot_t::SkipBlock(FILE *fp)
{
  int blocksize,blocksize2;
  
  ReadBlockSize(blocksize);	  
  fseek(fp, blocksize, SEEK_CUR);
  ReadBlockSize(blocksize2);
  assert(blocksize==blocksize2);
  return blocksize;
}
void Snapshot_t::GetFileName(Parameter_t &param, int ifile, string &filename)
{
  FILE *fp;
  char buf[1024];
  
  sprintf(buf,"%s/snapdir_%03d/%s_%03d.%d",param.SnapshotPath.c_str(),SnapshotId,param.SnapshotFileBase.c_str(),SnapshotId,ifile);
  if(ifile==0)
	if(!file_exist(buf)) sprintf(buf,"%s/%s_%03d",param.SnapshotPath.c_str(),param.SnapshotFileBase.c_str(),SnapshotId); //try the other convention
  if(!file_exist(buf)) sprintf(buf,"%s/%s_%03d.%d",param.SnapshotPath.c_str(),param.SnapshotFileBase.c_str(),SnapshotId,ifile); //try the other convention
  if(!file_exist(buf)) sprintf(buf,"%s/%d/%s.%d",param.SnapshotPath.c_str(),SnapshotId,param.SnapshotFileBase.c_str(),ifile);//for BJL's RAMSES output
  if(!file_exist(buf))
  {
	cerr<<"Failed to find a snapshot file at "<<buf<<endl;
	exit(1);
  }
  filename=buf;
}
bool Snapshot_t::ReadFileHeader(FILE *fp, SnapshotHeader_t &header)
{
  //read the header part, assign header extensions, and check byteorder
  int headersize(SNAPSHOT_HEADER_SIZE),headersize_byteswap(SNAPSHOT_HEADER_SIZE);
  swap_Nbyte(&headersize_byteswap,1,sizeof(headersize_byteswap));

  bool NeedByteSwap;
  int dummy,dummy2;
  fread(&dummy,sizeof(dummy),1,fp);
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
  myfread(&header.time,sizeof(double),1,fp);
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
	fprintf(stderr,"error!record brackets not match for header!\t%d,%d\n",dummy,dummy2);
	exit(1);
  } 
  
  return NeedByteSwap;
}
HBTInt Snapshot_t::ReadNumberOfDMParticles(Parameter_t & param, int ifile)
{
  FILE * fp;
  string filename;
  GetFileName(param, ifile, filename);
  myfopen(fp, filename.c_str(), "r");
  SnapshotHeader_t header;
  ReadFileHeader(fp, header);
  fclose(fp);
  return header.npart[1];
}
void Snapshot_t::LoadHeader(Parameter_t & param, int ifile)
{
  //read the header part, assign header extensions, and check byteorder

  FILE *fp; 
  string filename;
  GetFileName(param, ifile, filename);
  
  myfopen(fp, filename.c_str(), "r");
  NeedByteSwap=ReadFileHeader(fp, Header);
  
  if(Header.BoxSize!=param.BoxSize)
  {
	cerr<<"BoxSize not match input: read "<<Header.BoxSize<<"; expect "<<param.BoxSize<<endl;
	cerr<<"Maybe the length unit differ? Excpected unit: "<<param.LengthInMpch<<" Msol/h\n";
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
  
  if(param.SnapshotHasIdBlock)
  {
	blocksize=SkipIdBlock(fp);
	IntTypeSize=blocksize/NumPartInFile;
	assert(sizeof(int)==IntTypeSize||sizeof(long)==IntTypeSize);
	assert(sizeof(HBTInt)>=IntTypeSize);
  }
  
  fclose(fp);
  
  Header.Hz=PhysicalConst::H0 * sqrt(Header.Omega0 / (Header.time * Header.time * Header.time) 
  + (1 - Header.Omega0 - Header.OmegaLambda) / (Header.time * Header.time)
  + Header.OmegaLambda);//Hubble param for the current catalogue;
  
  //npartTotal is not reliable
  NumberOfParticles=0;
  for(int iFile=0;iFile<Header.num_files;iFile++)
  {
	int n=ReadNumberOfDMParticles(param, iFile);
	NumberOfDMParticleInFiles.push_back(n);
	OffsetOfDMParticleInFiles.push_back(NumberOfParticles);
	NumberOfParticles+=n;
  }
  
}

void Snapshot_t::Load(Parameter_t & param, int snapshot_index, bool load_id, bool load_position, bool load_velocity, bool load_mass, bool fill_particle_hash)
{ 
  SetSnapshotIndex(param, snapshot_index);
  PeriodicBox=param.PeriodicBoundaryOn;
  LoadHeader(param);
  if(load_id)
  {
	LoadId(param);
	if(fill_particle_hash)
	  FillParticleHash();
  }
  if(load_position)
	LoadPosition(param);
  if(load_velocity)
	LoadVelocity(param);
  if(load_mass)
	if(0.==Header.mass[1]) LoadMass(param);
}

size_t Snapshot_t::ReadBlock(FILE *fp, void *block, const size_t n_read, const size_t n_skip_before, const size_t n_skip_after)
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
void * Snapshot_t::LoadBlock(Parameter_t &param, int block_id, size_t element_size, int dimension, bool is_massblock)
{
  char * buf=static_cast<char *>(::operator new(element_size*NumberOfParticles*dimension));
  //#pragma omp parallel //can do task parallelization here.
  for(int iFile=0; iFile<Header.num_files; iFile++)
  {
	FILE *fp; 
	string filename;
	GetFileName(param, iFile, filename);
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

void Snapshot_t::LoadId(Parameter_t &param)
{//only do this after LoadHeader().
  
  if(!param.SnapshotHasIdBlock) return;
  
  void * buf=LoadBlock(param, BLOCK_ID_PID, IntTypeSize, 1);
  
  if(param.ParticleIdRankStyle)
  {
	cerr<<"Error: ParticleIdRankStyle not implemented yet\n";
	exit(1);
  }
  else
  {
	if(param.SnapshotIdUnsigned)//unsigned int
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

void Snapshot_t::LoadPosition(Parameter_t &param)
{
  void * buf=LoadBlock(param, BLOCK_ID_POS, RealTypeSize, 3);
  
  if(sizeof(HBTReal)==RealTypeSize)
	ComovingPosition=static_cast<HBTxyz *>(buf);
  else
  {
	ComovingPosition=static_cast<HBTxyz *>(::operator new(sizeof(HBTxyz)*NumberOfParticles));
	if(sizeof(HBTReal)==RealTypeSize)
	  AssignBlock(reinterpret_cast<HBTReal *>(ComovingPosition), static_cast<float *>(buf), NumberOfParticles*3);
	else
	  AssignBlock(reinterpret_cast<HBTReal *>(ComovingPosition), static_cast<double *>(buf), NumberOfParticles*3);
	::operator delete(buf);
  }
  
  if(param.PeriodicBoundaryOn)//regularize coord
  {
	for(HBTInt i=0;i<NumberOfParticles;i++)
	  for(int j=0;j<3;j++)
		ComovingPosition[i][j]=position_modulus(ComovingPosition[i][j], param.BoxSize);
  }
}

void Snapshot_t::LoadVelocity(Parameter_t &param)
{
  void * buf=LoadBlock(param, BLOCK_ID_VEL, RealTypeSize, 3);
  
  if(sizeof(HBTReal)==RealTypeSize)
	PhysicalVelocity=static_cast<HBTxyz *>(buf);
  else
  {
	PhysicalVelocity=static_cast<HBTxyz *>(::operator new(sizeof(HBTxyz)*NumberOfParticles));
	if(sizeof(HBTReal)==RealTypeSize)
	  AssignBlock(reinterpret_cast<HBTReal *>(PhysicalVelocity), static_cast<float *>(buf), NumberOfParticles*3);
	else
	  AssignBlock(reinterpret_cast<HBTReal *>(PhysicalVelocity), static_cast<double *>(buf), NumberOfParticles*3);
	::operator delete(buf);
  }
  
  HBTReal velocity_scale=sqrt(Header.time);
  for(HBTInt i=0;i<NumberOfParticles;i++)
	for(int j=0;j<3;j++)
	  PhysicalVelocity[i][j]*=velocity_scale;
}

void Snapshot_t::LoadMass(Parameter_t &param)
{
  void * buf=LoadBlock(param, BLOCK_ID_MASS, RealTypeSize, 1, true);
  
  if(sizeof(HBTReal)==RealTypeSize)
	ParticleMass=static_cast<HBTReal *>(buf);
  else
  {
	ParticleMass=static_cast<HBTReal *>(::operator new(sizeof(HBTReal)*NumberOfParticles));
	if(sizeof(HBTReal)==RealTypeSize)
	  AssignBlock(ParticleMass, static_cast<float *>(buf), NumberOfParticles);
	else
	  AssignBlock(ParticleMass, static_cast<double *>(buf), NumberOfParticles);
	::operator delete(buf);
  }  
}

void Snapshot_t::Clear()
/*reset to empty*/
{
#define RESET(x) {::operator delete(x);x=nullptr;}
  RESET(ParticleId);
  RESET(ComovingPosition);
  RESET(PhysicalVelocity);
  RESET(ParticleMass);
#undef RESET 
  ParticleHash.clear();//even if you don't do this, the destructor will still clean up the memory.
//   cout<<NumberOfParticles<<" particles cleared from snapshot "<<SnapshotIndex<<endl;
  NumberOfParticles=0;
}

void Snapshot_t::FillParticleHash()
{
  cout<<"Filling Hash Table...\n";
  ParticleHash.rehash(NumberOfParticles);
  ParticleHash.reserve(NumberOfParticles);
  for(HBTInt i=0;i<NumberOfParticles;i++)
	ParticleHash[ParticleId[i]]=i;
//   cout<<ParticleHash.bucket_count()<<" buckets used; load factor "<<ParticleHash.load_factor()<<endl;
}
void Snapshot_t::ClearParticleHash()
{
  ParticleHash.clear();
}
void Snapshot_t::AveragePosition(HBTxyz& CoM, const IndexList_t & Particles) const
/*mass weighted average position*/
{
	ParticleIndex_t i,j,np=Particles.size();
	double sx[3],origin[3],msum;
	static bool Periodic=PeriodicBox;
	static HBTReal boxsize=Header.BoxSize, boxhalf=boxsize/2.;
	
	if(0==np) return;
	
	sx[0]=sx[1]=sx[2]=0.;
	msum=0.;
	if(Periodic)
	  for(j=0;j<3;j++)
		origin[j]=GetComovingPosition(Particles[0])[j];
	
	for(i=0;i<np;i++)
	{
	  HBTReal m=GetParticleMass(Particles[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
	  if(Periodic)
		  sx[j]+=NEAREST(GetComovingPosition(Particles[i])[j]-origin[j], boxhalf, boxsize)*m;
	  else
		  sx[j]+=GetComovingPosition(Particles[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	{
		sx[j]/=msum;
		if(Periodic) sx[j]+=origin[j];
		CoM[j]=sx[j];
	}
}

void Snapshot_t::AverageVelocity(HBTxyz& CoV, const Snapshot_t::IndexList_t& Particles) const
/*mass weighted average velocity*/
{
	ParticleIndex_t i,j,np=Particles.size();
	double sv[3],msum;
	
	if(0==np) return;
	
	sv[0]=sv[1]=sv[2]=0.;
	msum=0.;
	
	for(i=0;i<np;i++)
	{
	  HBTReal m=GetParticleMass(Particles[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
		sv[j]+=GetPhysicalVelocity(Particles[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	  CoV[j]=sv[j]/msum;
}



#ifdef TEST_SIMU_IO
#include "../config_parser.h"

int main(int argc, char **argv)
{
  HBTConfig.ParseConfigFile(argv[1]);
  Snapshot_t snapshot;
  snapshot.Load(HBTConfig, 0, true, true, false, true);
  cout<<snapshot.GetNumberOfParticles()<<endl;
  cout<<snapshot.GetParticleId(10)<<endl;
  cout<<snapshot.GetComovingPosition(10)<<endl;
  cout<<snapshot.GetParticleMass(10)<<','<<snapshot.GetParticleMass(100)<<endl;
  cout<<snapshot.GetParticleIndex(snapshot.GetParticleId(10))<<endl;
  return 0;
}
#endif
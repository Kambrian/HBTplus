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
  
  myfread(header.npart,sizeof(int),TypeMax,fp);
  myfread(header.mass,sizeof(double),TypeMax,fp);
  myfread(&header.ScaleFactor,sizeof(double),1,fp);
  myfread(&header.redshift,sizeof(double),1,fp);
  myfread(&header.flag_sfr,sizeof(int),1,fp);
  myfread(&header.flag_feedback,sizeof(int),1,fp);
  myfread(header.npartTotal,sizeof(int),TypeMax,fp);
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
HBTInt ParticleSnapshot_t::ReadNumberOfParticles(int ifile)
{
  FILE * fp;
  string filename;
  GetFileName( ifile, filename);
  myfopen(fp, filename.c_str(), "r");
  SnapshotHeader_t header;
  ReadFileHeader(fp, header);
  fclose(fp);
  HBTInt np=0;
  for(int i=0;i<TypeMax;i++)
	np+=header.npart[i];
  return np;
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
  for(int i=0;i<TypeMax;i++) 
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
	int n=ReadNumberOfParticles( iFile);
	NumberOfParticleInFiles.push_back(n);
	OffsetOfParticleInFiles.push_back(NumberOfParticles);
	NumberOfParticles+=n;
  }
  
}

void ParticleSnapshot_t::Load(int snapshot_index, bool fill_particle_hash)
{ 
  SetSnapshotIndex(snapshot_index);
  
  Particles.clear();
  LoadHeader();
  {
  HBTInt np=0;
  np=accumulate(NumberOfParticleInFiles.begin(), NumberOfParticleInFiles.end(), np);
  Particles.reserve(np);
  }
  
 for(int iFile=0; iFile<Header.num_files; iFile++)
	ReadFile(iFile);

  if(HBTConfig.PeriodicBoundaryOn)//regularize coord
  {
	for(HBTInt i=0;i<size();i++)
	  for(int j=0;j<3;j++)
		Particles[i].ComovingPosition[j]=position_modulus(Particles[i].ComovingPosition[j], HBTConfig.BoxSize);
  }
  
  HBTReal velocity_scale=sqrt(Header.ScaleFactor);
  for(HBTInt i=0;i<size();i++)
	for(int j=0;j<3;j++)
	  Particles[i].PhysicalVelocity[j]*=velocity_scale;
	
	cout<<" ( "<<Header.num_files<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
	
	if(fill_particle_hash)
	FillParticleHash();
}

#define ReadScalarBlock(dtype, Attr) {\
FortranBlock <dtype> block(fp, n_read, n_skip, NeedByteSwap);\
  for(HBTInt i=0;i<n_read;i++)\
	NewParticles[i].Attr=block[i];	\
}

#define ReadXYZBlock(dtype, Attr) {\
  FortranBlock <dtype> block(fp, n_read*3, n_skip*3, NeedByteSwap);\
  auto pblock=block.data_reshape();\
  for(HBTInt i=0;i<n_read;i++)\
	for(int j=0;j<3;j++)\
	  NewParticles[i].Attr[j]=pblock[i][j];\
}


#define ReadMassBlock(dtype) {\
size_t n_read_mass=0;\
vector <HBTInt> offset_mass(TypeMax);\
for(int itype=0;itype<TypeMax;itype++)\
	  if(MassDataPresent(itype))\
	  {\
		offset_mass[itype]=n_read_mass;\
		n_read_mass+=header.npart[itype];\
	  }\
FortranBlock <dtype> block(fp, n_read_mass, n_skip, NeedByteSwap);\
  for(int itype=0;itype<TypeMax;itype++)\
  {\
	auto p=NewParticles+offset[itype];\
	dtype *pblock=block.data()+offset_mass[itype];\
	if(MassDataPresent(itype))\
	  for(HBTInt i=0;i<header.npart[itype];i++)\
		p[i].Mass=pblock[i];\
	else\
	{\
	  auto m=header.mass[itype];\
	  for(HBTInt i=0;i<header.npart[itype];i++)\
		p[i].Mass=m;\
	}\
  }\
}

#define ReadEnergyBlock(dtype) {\
FortranBlock <dtype> block(fp, header.npart[0], 0, NeedByteSwap);\
  for(HBTInt i=0;i<n_read;i++)\
	NewParticles[i].InternalEnergy=block[i];	\
}

void ParticleSnapshot_t::ReadFile(int iFile)
{ 
  FILE *fp; 
  string filename;
  GetFileName( iFile, filename);
  myfopen(fp, filename.c_str(), "r");
  SnapshotHeader_t header;
  ReadFileHeader(fp, header);
  size_t n_read=accumulate(begin(header.npart), end(header.npart), (size_t)0), n_skip=0;
  vector <HBTInt> offset(TypeMax);
  CompileOffsets(begin(header.npart), end(header.npart), offset.begin());
  
  Particles.resize(Particles.size()+n_read);
  const auto NewParticles=Particles.end()-n_read;
  
	if(RealTypeSize==4)
	{
	  ReadXYZBlock(float, ComovingPosition)
	  ReadXYZBlock(float, PhysicalVelocity)
	}
	else
	{
	  ReadXYZBlock(double, ComovingPosition)
	  ReadXYZBlock(double, PhysicalVelocity)
	}

	if(HBTConfig.SnapshotHasIdBlock)
	{
	  if(HBTConfig.ParticleIdRankStyle)
	  {
		cout<<"Error: ParticleIdRankStyle not implemented yet\n";
		exit(1);
	  }
	
	  if(IntTypeSize==4)
	  {
		if(HBTConfig.SnapshotIdUnsigned)//unsigned int
		  ReadScalarBlock(unsigned, Id)
		else
		  ReadScalarBlock(int, Id)
	  }
	  else
		ReadScalarBlock(long, Id)
	}
	else
	{
	  HBTInt id_now=OffsetOfParticleInFiles[iFile];
	  for(HBTInt i=0;i<n_read;i++)
		NewParticles[i].Id=id_now+i;
	}
  
  #define MassDataPresent(i) ((0==header.mass[i])&&(header.npartTotal[i]))
  if(RealTypeSize==4)
	ReadMassBlock(float)
  else
	ReadMassBlock(double)
#undef MassDataPresent
  
  if(RealTypeSize==4)
	ReadEnergyBlock(float)
  else
	ReadEnergyBlock(double)
  
  for(int itype=0;itype<TypeMax;++itype)
  {
	auto p=NewParticles+offset[itype];
	for(HBTInt i=0;i<header.npart[itype];i++)
	  p[i].Type=static_cast<ParticleType_t>(itype);
  }
	
  if(feof(fp))
  {
	cout<<"Error: End-of-File when reading "<<filename<<endl;
	exit(1);
  }
  
  fclose(fp);
}

void ParticleSnapshot_t::Clear()
/*reset to empty*/
{
  vector<Particle_t>().swap(Particles);
  ClearParticleHash();//even if you don't do this, the destructor will still clean up the memory.
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
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
#include "gadget_io.h"

void GadgetHeader_t::create_MPI_type(MPI_Datatype& dtype)
{
  /*to create the struct data type for communication*/	
  GadgetHeader_t &p=*this;
  #define NumAttr 13
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint   offsets[NumAttr], origin,extent;
  
  MPI_Get_address(&p,&origin);
  MPI_Get_address((&p)+1,&extent);//to get the extent of s
  extent-=origin;
  
  int i=0;
  #define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
  RegisterAttr(npart, MPI_INT, TypeMax)
  RegisterAttr(mass, MPI_DOUBLE, TypeMax)
  RegisterAttr(ScaleFactor, MPI_DOUBLE, 1)
  RegisterAttr(redshift, MPI_DOUBLE, 1)
  RegisterAttr(flag_sfr, MPI_INT, 1)
  RegisterAttr(flag_feedback, MPI_INT, 1)
  RegisterAttr(npartTotal, MPI_UNSIGNED, 1)
  RegisterAttr(flag_cooling, MPI_INT, 1)
  RegisterAttr(num_files, MPI_INT, 1)
  RegisterAttr(BoxSize, MPI_DOUBLE, 1)
  RegisterAttr(OmegaM0, MPI_DOUBLE, 1)
  RegisterAttr(OmegaLambda0, MPI_DOUBLE, 1)
  RegisterAttr(HubbleParam, MPI_DOUBLE, 1)
  #undef RegisterAttr
  assert(i==NumAttr);
  
  MPI_Type_create_struct(NumAttr,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef NumAttr
}

GadgetReader_t::GadgetReader_t(MpiWorker_t &world, int snapshot_id, vector< Particle_t >& particles, Cosmology_t& cosmology): SnapshotId(snapshot_id), Particles(particles), Cosmology(cosmology), Header()
{
  NeedByteSwap=false;
  IntTypeSize=0;
  RealTypeSize=0;
  Load(world);
}

#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
#define SkipPositionBlock(fp) SkipFortranBlock(fp, NeedByteSwap)
#define SkipVelocityBlock(fp) SkipFortranBlock(fp, NeedByteSwap)	 
#define SkipIdBlock(fp) SkipFortranBlock(fp, NeedByteSwap)
#define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)
void GadgetReader_t::GetGadgetFileName(int ifile, string &filename)
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
bool GadgetReader_t::ReadGadgetFileHeader(FILE *fp, GadgetHeader_t &header)
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

HBTInt GadgetReader_t::ReadGadgetNumberOfParticles(int ifile)
{
  FILE * fp;
  string filename;
  GetGadgetFileName( ifile, filename);
  myfopen(fp, filename.c_str(), "r");
  GadgetHeader_t header;
  ReadGadgetFileHeader(fp, header);
  fclose(fp);
  HBTInt np=0;
#ifdef DM_ONLY
  np=header.npart[TypeDM];
#else
  for(int i=0;i<TypeMax;i++)
	np+=header.npart[i];
#endif
  return np;
}

void GadgetReader_t::LoadGadgetHeader(int ifile)
{
  //read the header part, assign header extensions, and check byteorder

  FILE *fp; 
  string filename;
  GetGadgetFileName( ifile, filename);
  
  myfopen(fp, filename.c_str(), "r");
  NeedByteSwap=ReadGadgetFileHeader(fp, Header);
  
  if((HBTReal)Header.BoxSize!=HBTConfig.BoxSize)
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
// 	if(sizeof(HBTInt)<IntTypeSize)
// 	  cerr<<"WARNING: loading size "<<IntTypeSize<<" integer in snapshot with size "<<sizeof(HBTInt)<<" int in HBT. possible data overflow.\n Please use ./HBTdouble unless you know what you are doing.";
  }
  
  fclose(fp);
  
  //npartTotal is not reliable
  HBTInt np=0;
  for(int iFile=0;iFile<Header.num_files;iFile++)
  {
	int n=ReadGadgetNumberOfParticles( iFile);
	NumberOfParticleInFiles.push_back(n);
	OffsetOfParticleInFiles.push_back(np);
	np+=n;
  }
}

void GadgetReader_t::Load(MpiWorker_t &world)
{ 
  const int root=0;
  if(world.rank()==root)
    LoadGadgetHeader();
  MPI_Datatype MPI_GadgetHeader_t;
  GadgetHeader_t().create_MPI_type(MPI_GadgetHeader_t);
  MPI_Bcast(&Header,1, MPI_GadgetHeader_t, root, world.Communicator);
  MPI_Type_free(&MPI_GadgetHeader_t);
  world.SyncAtomBool(NeedByteSwap, root);
  world.SyncAtom(IntTypeSize, MPI_INT, root);
  world.SyncAtom(RealTypeSize, MPI_INT, root);
  world.SyncContainer(NumberOfParticleInFiles, MPI_HBT_INT, root);
  world.SyncContainer(OffsetOfParticleInFiles, MPI_HBT_INT, root);

  Cosmology.Set(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);  
#ifdef DM_ONLY
//   Cosmology.ParticleMass=Header.mass[TypeDM];
#endif
  
  int nfiles_skip, nfiles_end;
  AssignTasks(world.rank(), world.size(), Header.num_files, nfiles_skip, nfiles_end);
  {
  HBTInt np=0;
  np=accumulate(NumberOfParticleInFiles.begin()+nfiles_skip, NumberOfParticleInFiles.begin()+nfiles_end, np);
  Particles.reserve(np);
  }
  
  for(int i=0, ireader=0;i<world.size();i++, ireader++)
  {
	if(ireader==HBTConfig.MaxConcurrentIO) 
	{
	  ireader=0;//reset reader count
	  MPI_Barrier(world.Communicator);//wait for every thread to arrive.
	}
	if(i==world.rank())//read
	{
	  for(int iFile=nfiles_skip; iFile<nfiles_end; iFile++)
		ReadGadgetFile(iFile);
	}
  }
  
//   cout<<" ( "<<Header.num_files<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
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

#ifndef DM_ONLY
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
	  const dtype * pblock=block.data()+offset_mass[itype];\
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
#else
  #define ReadMassBlock(dtype) {\
  if(MassDataPresent(TypeDM))\
  {\
    size_t offset_mass=0;\
    for(int itype=0;itype<TypeDM;itype++)\
      if(MassDataPresent(itype)) offset_mass+=header.npart[itype];\
    FortranBlock <dtype> block(fp, n_read, offset_mass, NeedByteSwap);\
    for(HBTInt i=0;i<n_read;i++) NewParticles[i].Mass=block[i];\
  }\
  else\
  {\
    auto m=header.mass[TypeDM];\
    for(HBTInt i=0;i<n_read;i++) NewParticles[i].Mass=m;\
  }\
  }
#endif

#define ReadEnergyBlock(dtype) {\
FortranBlock <dtype> block(fp, header.npart[0], 0, NeedByteSwap);\
  for(HBTInt i=0;i<header.npart[0];i++)\
	NewParticles[i].InternalEnergy=block[i];	\
}

void GadgetReader_t::ReadGadgetFile(int iFile)
{ 
  FILE *fp; 
  string filename;
  GetGadgetFileName( iFile, filename);
  myfopen(fp, filename.c_str(), "r");
  GadgetHeader_t header;
  ReadGadgetFileHeader(fp, header);
#ifdef DM_ONLY
  size_t n_read=header.npart[TypeDM], n_skip=accumulate(header.npart, header.npart+TypeDM, size_t(0));
#else
  size_t n_read=accumulate(begin(header.npart), end(header.npart), (size_t)0), n_skip=0;
#endif
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

	if(HBTConfig.PeriodicBoundaryOn)//regularize coord
	{
	  HBTReal boxsize=HBTConfig.BoxSize;
	  for(HBTInt i=0;i<n_read;i++)
		for(int j=0;j<3;j++)
		  NewParticles[i].ComovingPosition[j]=position_modulus(NewParticles[i].ComovingPosition[j], boxsize);
	}
  
	HBTReal velocity_scale=sqrt(Header.ScaleFactor);
	for(HBTInt i=0;i<n_read;i++)
	  for(int j=0;j<3;j++)
		NewParticles[i].PhysicalVelocity[j]*=velocity_scale;
	
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

#ifndef DM_ONLY
#ifdef HAS_THERMAL_ENERGY
  if(RealTypeSize==4)
	ReadEnergyBlock(float)
  else
	ReadEnergyBlock(double)
#endif
	
  for(int itype=0;itype<TypeMax;++itype)
  {
	auto p=NewParticles+offset[itype];
	for(HBTInt i=0;i<header.npart[itype];i++)
	  p[i].Type=static_cast<ParticleType_t>(itype);
  }
#endif

  if(feof(fp))
  {
	cout<<"Error: End-of-File when reading "<<filename<<endl;
	exit(1);
  }
  
  fclose(fp);
}
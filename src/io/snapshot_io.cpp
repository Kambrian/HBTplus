using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <chrono>
#include <numeric> 
#include <cstdlib>
#include <cstdio>
#include <boost/concept_check.hpp>

#include "../mpi_wrapper.h"
#include "../snapshot.h"
#include "../mymath.h"

#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
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
	cerr<<"Error: record brackets not match for header!"<<dummy<<dummy2<<endl;
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
  
  if((HBTReal)Header.BoxSize!=HBTConfig.BoxSize)
  {
	cerr<<"BoxSize not match input: read "<<Header.BoxSize<<"; expect "<<HBTConfig.BoxSize<<endl;
	cerr<<"Maybe the length unit differ? Excpected unit: "<<HBTConfig.LengthInMpch<<" Msol/h\n";
	exit(1);
  }
  
  //set datatypes
  HBTInt NumPartInFile=0;
  for(int i=0;i<NUMBER_OF_PARTICLE_TYPES;i++) 
	NumPartInFile+=Header.npart[i];
  int blocksize=SkipBlock(fp);
  RealTypeSize=blocksize/NumPartInFile/3;
  assert(sizeof(float)==RealTypeSize||sizeof(double)==RealTypeSize);
  blocksize=SkipBlock(fp);
  assert(blocksize==RealTypeSize*NumPartInFile*3);
  if(sizeof(HBTReal)<RealTypeSize)
	cerr<<"WARNING: loading size "<<RealTypeSize<<" float in snapshot with size "<<sizeof(HBTReal)<<" float in HBT. possible loss of accuracy.\n Please use ./HBTdouble unless you know what you are doing.";
  
  if(HBTConfig.SnapshotHasIdBlock)
  {
	blocksize=SkipBlock(fp);
	IntTypeSize=blocksize/NumPartInFile;
	assert(sizeof(int)==IntTypeSize||sizeof(long)==IntTypeSize);
	assert(sizeof(HBTInt)>=IntTypeSize);
	if(sizeof(HBTInt)<IntTypeSize)
	  cerr<<"WARNING: loading size "<<IntTypeSize<<" integer in snapshot with size "<<sizeof(HBTInt)<<" int in HBT. possible data overflow.\n Please use ./HBTdouble unless you know what you are doing.";
  }
  
  fclose(fp);
  
  SetEpoch(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);
  
  //npartTotal is not reliable
  HBTInt np_allfiles=0;
  for(int iFile=0;iFile<Header.num_files;iFile++)
  {
	int n=ReadNumberOfDMParticles( iFile);
	NumberOfDMParticleInFiles.push_back(n);
	OffsetOfDMParticleInFiles.push_back(np_allfiles);
	np_allfiles+=n;
  }
  NumberOfParticlesOnAllNodes=np_allfiles;
}

void ParticleSnapshot_t::Load(MpiWorker_t & world, int snapshot_index, bool fill_particle_hash)
{ 
  Clear();
  SetSnapshotIndex(snapshot_index);
   
  {//load header
  const int root=0;
  if(world.rank()==root)
	LoadHeader(0);
  MPI_Datatype MPI_SnapshotHeader_t;
  SnapshotHeader_t().create_MPI_type(MPI_SnapshotHeader_t);
  MPI_Bcast(&Header,1, MPI_SnapshotHeader_t, root, world.Communicator);
  MPI_Type_free(&MPI_SnapshotHeader_t);
  world.SyncAtom(OmegaM0, MPI_HBT_REAL, root);
  world.SyncAtom(OmegaLambda0, MPI_HBT_REAL, root);
  world.SyncAtom(Hz, MPI_HBT_REAL, root);
  world.SyncAtom(ScaleFactor, MPI_HBT_REAL, root);
  world.SyncAtomBool(NeedByteSwap, root);
  world.SyncAtom(IntTypeSize, MPI_INT, root);
  world.SyncAtom(RealTypeSize, MPI_INT, root);
  world.SyncContainer(NumberOfDMParticleInFiles, MPI_HBT_INT, root);
  world.SyncContainer(OffsetOfDMParticleInFiles, MPI_HBT_INT, root);
  world.SyncAtom(NumberOfParticlesOnAllNodes, MPI_HBT_INT, root);
  }
  
  int nfiles_skip, nfiles_end;
  AssignTasks(world.rank(), world.size(), Header.num_files, nfiles_skip, nfiles_end);
  {
  HBTInt np=0;
  np=accumulate(NumberOfDMParticleInFiles.begin()+nfiles_skip, NumberOfDMParticleInFiles.begin()+nfiles_end, np);
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
		ReadFile(iFile);
	}
  }
      
  if(!HBTConfig.SnapshotHasIdBlock)
  for(HBTInt i=0;i<size();i++)
	Particles[i].Id=OffsetOfDMParticleInFiles[nfiles_skip]+i;
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
	
	/*
	 i f((size()>1)&&(ParticleMass[0]!=ParticleMass[1]))
	 {
	 cout<<"Error: DM particles have different mass? not supported!\n";
	 exit(1);
}
else
  Header.mass[1]=ParticleMass[0];*/
	
	cout<<"Finished Reading on thread "<<world.rank()<<" from file "<<nfiles_skip<<" to "<<nfiles_end-1;
	cout<<" ( "<<Header.num_files<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
	
	ExchangeParticles(world);
	
	if(fill_particle_hash)
	FillParticleHash();
	
	cout<<"IdRange=("<<IdMin<<","<<IdMax<<") on "<<world.rank()<<endl;
}
inline int GetGrid(HBTReal x, HBTReal step, int dim)
{
  int i=floor(x/step);
  if(i<0) i=0;
  if(i>=dim) i=dim-1;
  return i;
}
inline int AssignCell(HBTxyz & Pos, const HBTReal step[3], const vector <int> &dims)
{
  #define GRIDtoRank(g0,g1,g2) (((g0)*dims[1]+(g1))*dims[2]+(g2))
  #define GID(i) GetGrid(Pos[i], step[i], dims[i])
  return GRIDtoRank(GID(0), GID(1), GID(2));
}
/* chunked version
void ParticleSnapshot_t::ExchangeParticles(MpiWorker_t &world)
{
  #define MSG_SIZE 1024
  
  auto dims=ClosestFactors(world.size(), 3);
  HBTReal step[3];
  for(int i=0;i<3;i++)
	step[i]=HBTConfig.BoxSize/dims[i];
  
  typedef vector <Particle_t> ParticleList_t;
  ParticleList_t NewParticles;
  vector <ParticleList_t> SendCells(world.size());
  vector <ParticleList_t> ReceiveCells(world.size());
  auto current_particle=Particles.begin(), end_particle=current_particle, max_particle=Particles.end();
  HBTInt nloop=ceil(1.*Particles.size()/MSG_SIZE);
  nloop=mpi::all_reduce(world, nloop, mpi::maximum<HBTInt>());
  for(HBTInt iloop=0;iloop<nloop;iloop++)
  {
	end_particle=min(current_particle+MSG_SIZE, max_particle);
	for(;current_particle<end_particle;current_particle++)//pick
	{
	  int rank=AssignCell(current_particle->ComovingPosition, step, dims);
	  SendCells[rank].push_back(*current_particle);
	}
	all_to_all(world, SendCells, ReceiveCells);//deliver
  	for(int i=0;i<world.size();i++)//insert
	{
	  NewParticles.insert(NewParticles.end(), ReceiveCells[i].begin(), ReceiveCells[i].end());
	  ReceiveCells[i].clear();
	  SendCells[i].clear();
	}
  }
  cout<<NewParticles.size()<<" particles received on node "<<world.rank()<<endl;
  Particles.swap(NewParticles);
} */
void ParallelStride(MpiWorker_t &world, vector <Particle_t> &Particles, HBTInt &offset, HBTInt steps)
{  
  while(steps)
  {
	HBTInt nmax=Particles.size()-offset;
	MPI_Comm newcomm;
	MPI_Comm_split(world.Communicator, nmax==0, 0, &newcomm);
	int newcomm_size;
	MPI_Comm_size(newcomm, &newcomm_size);
	
	HBTInt n=0;
	if(nmax)
	{
	  n=steps/newcomm_size;
	  if(n*newcomm_size<steps) n++;
	  if(n>nmax) n=nmax;

	  HBTInt pid=Particles[offset+n-1].Id;
	  HBTInt MinId;
	  MPI_Allreduce(&pid, &MinId, 1, MPI_HBT_INT, MPI_MIN, newcomm);
	
	  if(pid>MinId)//fastforward sparse ranks
	  {
		n=offset;
		while(Particles[offset].Id<MinId)
		  offset++;
		n=offset-n;
	  }
	  else
		offset+=n;
	}
	
	MPI_Comm_free(&newcomm);
	
	MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
	steps-=n;
  }
}
bool ParticleSnapshot_t::IsContiguousId(MpiWorker_t &world, HBTInt &GlobalIdMin)
{
  MPI_Comm newcomm;
  MPI_Comm_split(world.Communicator, Particles.size()==0, 0, &newcomm);
  int newcomm_size, newrank;
  MPI_Comm_size(newcomm, &newcomm_size);
  MPI_Comm_rank(newcomm, &newrank);
  
  int flag_contig=0;
  if(Particles.size())
  {
	IdMin=Particles.front().Id;
	IdMax=Particles.back().Id;
	MPI_Reduce(MPI_IN_PLACE, &IdMin, 1, MPI_HBT_INT, MPI_MIN, 0, newcomm);
	MPI_Reduce(MPI_IN_PLACE, &IdMax, 1, MPI_HBT_INT, MPI_MAX, 0, newcomm);
	if(newrank) 
	  IdMin=0;
	else
	{
	  flag_contig=(IdMax-IdMin+1==NumberOfParticlesOnAllNodes);
	  if(flag_contig) cout<<"Contiguous particle Ids."<<endl;
	}
  }
  else
	IdMin=0;
  
  MPI_Comm_free(&newcomm);
  MPI_Allreduce(MPI_IN_PLACE, &flag_contig, 1, MPI_INT, MPI_LOR, world.Communicator);
  MPI_Allreduce(&IdMin, &GlobalIdMin, 1, MPI_HBT_INT, MPI_BOR, world.Communicator);
  
  return flag_contig;
}

void ParticleSnapshot_t::PartitionParticles(MpiWorker_t &world, vector <int> &offset)
{
  int nremainder=NumberOfParticlesOnAllNodes%world.size();
  HBTInt nnew=NumberOfParticlesOnAllNodes/world.size()+1;
  
  HBTInt GlobalIdMin;
  if(IsContiguousId(world, GlobalIdMin))
  {
	HBTInt & upperbound=GlobalIdMin;
	int rank=0, pid=0;
	while(pid<Particles.size())
	{
	  if(Particles[pid].Id<upperbound) 
		pid++;
	  else
	  {
		offset[rank]=pid;
		if(rank==nremainder) 
		  nnew--;
		rank++;
		upperbound+=nnew;
	  }
	}
	while(rank<world.size())
	  offset[rank++]=pid;
	assert(pid==Particles.size());
  }
  else 
  {
	HBTInt this_offset=0;
	for(int i=0;i<world.size();i++)
	{
	  offset[i]=this_offset;
	  if(i==nremainder) 
		nnew--;
	  ParallelStride(world, Particles, this_offset, nnew);
	}
	assert(this_offset==Particles.size());
  }
}
inline bool CompParticleId(const Particle_t &a, const Particle_t &b)
{
  return a.Id<b.Id;
}
void ParticleSnapshot_t::ExchangeParticles(MpiWorker_t &world)
{
  sort(Particles.begin(), Particles.end(), CompParticleId);
  
  vector <int> SendOffsets(world.size()+1), SendSizes(world.size(), 0);
  PartitionParticles(world, SendOffsets);
  SendOffsets.back()=Particles.size();
  for(int i=0;i<world.size();i++)
	SendSizes[i]=SendOffsets[i+1]-SendOffsets[i];
  
  vector <int> ReceiveSizes(world.size(),0), ReceiveOffsets(world.size());
  MPI_Alltoall(SendSizes.data(), 1, MPI_INT, ReceiveSizes.data(), 1, MPI_INT, world.Communicator);
  vector <Particle_t> ReceivedParticles;
  ReceivedParticles.resize(CompileOffsets(ReceiveSizes, ReceiveOffsets));
  
  MPI_Datatype MPI_HBT_Particle;
  Particle_t().create_MPI_type(MPI_HBT_Particle);
  MPI_Alltoallv(Particles.data(), SendSizes.data(), SendOffsets.data(), MPI_HBT_Particle, 
				ReceivedParticles.data(), ReceiveSizes.data(), ReceiveOffsets.data(), MPI_HBT_Particle, world.Communicator);

  MPI_Type_free(&MPI_HBT_Particle);
  
  Particles.swap(ReceivedParticles);
  
  cout<<Particles.size()<<" particles received on node "<<world.rank()<<endl;
}

#define ReadScalarBlock(dtype, Attr) {\
FortranBlock <dtype> block(fp, n_read, n_skip, NeedByteSwap);\
  auto p=Particles.data()+Particles.size()-n_read;\
  for(HBTInt i=0;i<n_read;i++)\
	p[i].Attr=block[i];	\
}

#define ReadXYZBlock(dtype, Attr) {\
  FortranBlock <dtype> block(fp, n_read*3, n_skip*3, NeedByteSwap);\
  auto p=Particles.data()+Particles.size()-n_read;\
  auto pblock=block.data_reshape();\
  for(HBTInt i=0;i<n_read;i++)\
	for(int j=0;j<3;j++)\
	  p[i].Attr[j]=pblock[i][j];\
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
  
  Particles.resize(Particles.size()+n_read);
  
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
  
  #define MassDataPresent(i) ((0==header.mass[i])&&(header.npartTotal[i]))
  if(MassDataPresent(1))
  {
	size_t n_skip_old=n_skip;
	if(!MassDataPresent(0)) n_skip=0;
	
	if(RealTypeSize==4)
	  ReadScalarBlock(float, Mass)
	else
	  ReadScalarBlock(double, Mass)
	
	n_skip=n_skip_old;//restore 
  }
  else
  {
	for(auto it=Particles.end()-n_read; it<Particles.end();it++)
	  it->Mass=header.mass[1]; 

// 	if(!HBTConfig.SnapshotNoMassBlock)
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

void ParticleSnapshot_t::Clear()
/*reset to empty*/
{
  #define RESET(x, T) {vector <T>().swap(x);}
  RESET(Particles, Particle_t);
  #undef RESET 
  ClearParticleHash();//even if you don't do this, the destructor will still clean up the memory.
  //   cout<<NumberOfParticles<<" particles cleared from snapshot "<<SnapshotIndex<<endl;
  NumberOfParticlesOnAllNodes=0;
}

#ifdef TEST_snapshot_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  mpi::environment env;
  MpiWorker_t world;
  #ifdef _OPENMP
  omp_set_nested(0);
  #endif
  
  int snapshot_start, snapshot_end;
  if(0==world.rank())
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  HBTConfig.BroadCast(world, 0, snapshot_start, snapshot_end);
  
  ParticleSnapshot_t snapshot;
  snapshot.Load(world, snapshot_start, true);
  
  cout<<snapshot.size()<<" particles loaded on thread "<<world.rank()<<endl;
  cout<<"Particle 10: "<<snapshot.GetId(10);
  cout<<snapshot.GetComovingPosition(10);
  cout<<snapshot.GetMass(10)<<','<<snapshot.GetMass(100)<<endl;
  cout<<snapshot.GetIndex(snapshot.GetId(10))<<endl<<endl;
  
  return 0;
}
#endif
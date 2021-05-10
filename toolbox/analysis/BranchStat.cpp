/*to count the number of particles in each direct infall branch, for Wenting's oPDF analysis*/
using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
//#include <boost/concept_check.hpp>

#include "../../src/hdf_wrapper.h"
#include "../../src/mpi_wrapper.h"
#include "../../src/subhalo.h"
#include "../../src/halo.h"
#include "../../src/mymath.h"

typedef long BranchIdType_t;
#define MPI_BranchIdType MPI_LONG
#define H5T_BranchId H5T_NATIVE_LONG 
#define BranchIDSnapBase 10000000000000000 //BranchId=HaloId+SnapshotId*BranchIDSnapBase

struct ParticleBranch_t
{
  HBTInt ParticleId;
  BranchIdType_t BranchId;
  int Index;//current order
  int Index0;//original order
  ParticleBranch_t(): ParticleId(0), BranchId(-1)
  {
  }
  void create_MPI_dtype(MPI_Datatype &dtype);
};
MPI_Datatype MPI_ParticleBranch_t;
void ParticleBranch_t::create_MPI_dtype(MPI_Datatype &MPI_dtype)
{
ParticleBranch_t &p=*this;
#define NumAttr 4
MPI_Datatype oldtypes[NumAttr];
int blockcounts[NumAttr];
MPI_Aint   offsets[NumAttr], origin,extent;

MPI_Get_address(&p,&origin);
MPI_Get_address((&p)+1,&extent);//to get the extent of s
extent-=origin;

int i=0;
#define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
RegisterAttr(ParticleId, MPI_HBT_INT, 1)
RegisterAttr(BranchId, MPI_BranchIdType, 1)
RegisterAttr(Index, MPI_INT, 1)
RegisterAttr(Index0, MPI_INT, 1)
#undef RegisterAttr
assert(i==NumAttr);

MPI_Type_create_struct(NumAttr,blockcounts,offsets,oldtypes, &MPI_dtype);//some padding is added automatically by MPI as well
MPI_Type_create_resized(MPI_dtype,(MPI_Aint)0, extent, &MPI_dtype);
MPI_Type_commit(&MPI_dtype);
#undef NumAttr
}
inline bool CompBranch(const ParticleBranch_t & a, const ParticleBranch_t & b)
{
  return (a.BranchId>b.BranchId);
}
inline bool CompIndex0(const ParticleBranch_t & a, const ParticleBranch_t & b)
{
  return (a.Index0<b.Index0);
}
struct TargetHalo_t
{
  vector <ParticleBranch_t> Particle;
  HBTInt Nfound;
  HBTInt NBranch;
  void Bcast(MpiWorker_t &world, int root)
  {
	int np=Particle.size();
	MPI_Bcast(&np, 1, MPI_INT, root, world.Communicator);
	if(world.rank()!=root)
	  Particle.resize(np);
	MPI_Bcast(Particle.data()+Nfound, Particle.size()-Nfound, MPI_ParticleBranch_t, root, world.Communicator);
  }
  void FindHost(HaloSnapshot_t &halosnap)
  {
	for(HBTInt ipart=Nfound;ipart<Particle.size();ipart++)
	{
	  auto &p=Particle[ipart];
	  p.BranchId=halosnap.ParticleHash.GetIndex(p.ParticleId);
      p.BranchId=halosnap.Halos[p.BranchId].HaloId;//bugfix: change to global HaloId instead of local
      if(p.BranchId!=SpecialConst::NullHaloId)
          p.BranchId+=halosnap.GetSnapshotId()*BranchIDSnapBase;//prepend snapshotId
	  p.Index=ipart;//new index
// 		BranchCount[p.BranchId]++;
	}
  }
  void ReduceHost(MpiWorker_t &world)
  {//master
	unordered_map<BranchIdType_t, HBTInt> BranchCounter;
	vector <bool> FlagWait(world.size(), true);
	vector <ParticleBranch_t> ParticleBuffer(Particle.size());
	int nwait=world.size()-1;
	while(nwait>0)
	{
	  for(int rank=1;rank<world.size();rank++)
	  {
		if(FlagWait[rank])
		{
		  int flagready=0;
		  MPI_Status stat;
		  MPI_Iprobe(rank, 0, world.Communicator, &flagready, &stat);
		  if(flagready)
		  {
			int count;
			MPI_Get_count(&stat, MPI_ParticleBranch_t, &count);
			MPI_Recv(ParticleBuffer.data(), count, MPI_ParticleBranch_t, rank, 0, world.Communicator, &stat);
			FlagWait[rank]=false;
			nwait--;
			
			for(int i=0;i<count;i++)
			{
			  auto &h=ParticleBuffer[i];
			  Particle[h.Index].BranchId=h.BranchId;
			  BranchCounter[h.BranchId]++;
			}
		  }
		}
	  }
	}
	
	  /*find main branch*/
	  HBTInt MaxBId=-1,MaxCount=0, NewFound=0;
	  for(auto &x: BranchCounter)
	  {
		if(x.first!=SpecialConst::NullHaloId)
		{
		  NewFound+=x.second;
		  if(x.second>MaxCount)
		  {
			MaxCount=x.second;
			MaxBId=x.first;
		  }
		}
	  }
	  NewFound-=MaxCount;
	  HBTInt NumNewBranches=BranchCounter.size()-1;//uncount background
	  if(MaxBId!=-1) NumNewBranches--;//uncount main progenitor
	  /*remove main branch*/
	  for(HBTInt ipart=Nfound;ipart<Particle.size();ipart++)
	  {
		auto &p=Particle[ipart];
		if(p.BranchId==MaxBId) p.BranchId=-1;
	  }
	  /*sort according to branchid, move unfound to the end*/
	  sort(Particle.begin()+Nfound, Particle.end(), CompBranch);
	  Nfound+=NewFound;
	  NBranch+=NumNewBranches;
  }
  void SendHost(MpiWorker_t &world)
  {
	sort(Particle.begin()+Nfound, Particle.end(), CompBranch);//move -1 to end
	int iend, count;
	for(iend=Particle.size()-1;iend>=Nfound;iend--)
	  if(Particle[iend].BranchId!=-1) break;
	count=iend-Nfound+1;//bugfix
	MPI_Request Req;
	MPI_Isend(Particle.data()+Nfound, count, MPI_ParticleBranch_t, 0, 0, world.Communicator, &Req);
  }
};
class HaloList_t
{
public:
  vector <TargetHalo_t> Halo;
  void LoadSingle(int ihalo);
  void Load();
  void Save();
  void SaveSingle(int ihalo);
  void Bcast(MpiWorker_t &world, int root)
  {
	int nhalo=Halo.size();
	MPI_Bcast(&nhalo, 1, MPI_INT, root, world.Communicator);
	if(world.rank()!=0)
	  Halo.resize(nhalo);
	for(auto &&h: Halo)
	  h.Bcast(world, root);
  }
};

int main(int argc, char **argv)
{
    assert(sizeof(long)==16);//make sure BranchId dtype is long enough
   MPI_Init(&argc, &argv);
 MpiWorker_t world(MPI_COMM_WORLD);
#ifdef _OPENMP
 omp_set_nested(0);
#endif
 MPI_Comm MasterSlaveComm;
 MPI_Comm_split( MPI_COMM_WORLD, world.rank() == 0, world.rank(), &MasterSlaveComm);
 MpiWorker_t worker(MasterSlaveComm);
 
 ParticleBranch_t().create_MPI_dtype(MPI_ParticleBranch_t);
 
  int snapshot_start, snapshot_end;
  if(0==world.rank())
  {
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
	mkdir(HBTConfig.SubhaloPath.c_str(), 0755);
	HBTConfig.DumpParameters();
  }
  HBTConfig.BroadCast(world, 0, snapshot_start, snapshot_end);
  
  HaloList_t TargetHalos;
  if(world.rank()==0)
	TargetHalos.Load();
  for(int isnap=snapshot_end;isnap>=snapshot_start;isnap--)
  {
	TargetHalos.Bcast(world, 0);
	if(world.rank()>0)//slaves
	{
	  HaloSnapshot_t halosnap;
	  halosnap.Load(worker, isnap);
	  halosnap.FillParticleHash();
	  for(auto && halo: TargetHalos.Halo)
	  {
		halo.FindHost(halosnap);
		halo.SendHost(world);
	  }
	}
	else
	{
	  cout<<"Snap"<<isnap;
	  for(int i=0;i<TargetHalos.Halo.size();i++)
	  {
		if(i%10==0) cout<<".";
		TargetHalos.Halo[i].ReduceHost(world);
	  }
	  cout<<endl;
	}
  }
  
  TargetHalos.Save();
  
  MPI_Type_free(&MPI_ParticleBranch_t);
  MPI_Finalize();
  
  return 0;
}

void HaloList_t::Load()
{
  const int nhalos=944;
  Halo.resize(nhalos);
  for(int i=0;i<Halo.size();i++)
	LoadSingle(i);
}
void HaloList_t::LoadSingle(int ihalo)
{
  stringstream filename;
  filename<<"/cosma/home/durham/jvbq85/data/HBT/data/Millennium2/subcat.nostrip/postproc/IsolatedMWHalos.hdf5";
  stringstream dsetname;
  dsetname<<"Halo"<<setw(3)<<setfill('0')<<ihalo<<"/ParticleId";
  cout<<"Reading "<<filename.str()<<"/"<<dsetname.str()<<endl;
  hid_t dset, file=H5Fopen(filename.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hsize_t dims[2];
  dset=H5Dopen2(file, dsetname.str().c_str(), H5P_DEFAULT);
  int ndim=GetDatasetDims(dset, dims);
  vector <HBTInt> Particles(dims[0]);
  H5Dread(dset, H5T_HBTInt, H5S_ALL, H5S_ALL, H5P_DEFAULT, Particles.data());
  H5Dclose(dset);
  H5Fclose(file);
  
  Halo[ihalo].Particle.resize(Particles.size());
  for(HBTInt i=0;i<Particles.size();i++)
  {
	Halo[ihalo].Particle[i].ParticleId=Particles[i];
	Halo[ihalo].Particle[i].Index=i;
    Halo[ihalo].Particle[i].Index0=i;
  }
}
void HaloList_t::Save()
{
  for(int i=0;i<Halo.size();i++)
	SaveSingle(i);
}
void HaloList_t::SaveSingle(int ihalo)
{
  stringstream filename;
  filename<<"/cosma/home/durham/jvbq85/data/HBT/data/Millennium2/subcat.nostrip/postproc/IsolatedMWHalos.hdf5";
  stringstream dsetname;
  dsetname<<"Halo"<<setw(3)<<setfill('0')<<ihalo;
  cout<<"Saving "<<filename.str()<<"/"<<dsetname.str()<<endl;
  hid_t dset, file, grp;
  file=H5Fopen(filename.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  grp=H5Gopen(file, dsetname.str().c_str(), H5P_DEFAULT);
  auto &Particle=Halo[ihalo].Particle;
  sort(Particle.begin(), Particle.end(), CompIndex0);
  hsize_t dims[1]={Particle.size()}, ndim=1;
//   vector <HBTInt> ParticleId(Particle.size());
  vector <BranchIdType_t> BranchId(Particle.size());
  for(HBTInt i=0;i<Particle.size();i++)
  {
// 	ParticleId[i]=Particle[i].ParticleId;
	BranchId[i]=Particle[i].BranchId;
  }
//   writeHDFmatrix(file, ParticleId.data(), "dmid", ndim, dims, H5T_HBTInt); 
  writeHDFmatrix(grp, BranchId.data(), "BranchId", ndim, dims, H5T_BranchId); 
  dims[0]=1;
  writeHDFmatrix(grp, &Halo[ihalo].NBranch, "NBranch", ndim, dims, H5T_HBTInt);
  writeHDFmatrix(grp, &Halo[ihalo].Nfound, "TotNumPartInBranch", ndim, dims, H5T_HBTInt);
  cout<<"Nbranch="<<Halo[ihalo].NBranch<<", Nfound="<<Halo[ihalo].Nfound<<endl;
  H5Gclose(grp);
  H5Fclose(file);
}

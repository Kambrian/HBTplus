#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <glob.h>
#include <climits>
#include <algorithm>
#include <chrono>

#include "mymath.h"
#include "halo.h"
#include <mpi.h>

void create_MPI_Halo_Id_type(MPI_Datatype &MPI_HBTHalo_Id_t)
{
/*to create the struct data type for communication*/	
Halo_t p;
#define NumAttr 1
MPI_Datatype oldtypes[NumAttr]={MPI_HBT_INT};
int blockcounts[NumAttr]={1};
MPI_Aint   offsets[NumAttr], origin,extent;

MPI_Get_address(&p,&origin);
MPI_Get_address(&p.HaloId,offsets);
// MPI_Get_address(p.ComovingPosition.data(),offsets+1);//caution: this might be implementation dependent??
// MPI_Get_address(p.PhysicalVelocity.data(),offsets+2);
MPI_Get_address((&p)+1,&extent);//to get the extent of s

for(int i=0;i<NumAttr;i++)
  offsets[i]-=origin;

extent-=origin;

// assert(offsets[2]-offsets[1]==sizeof(HBTReal)*3);//to make sure HBTxyz is stored locally.

MPI_Type_create_struct(NumAttr,blockcounts,offsets,oldtypes, &MPI_HBTHalo_Id_t);//some padding is added automatically by MPI as well
MPI_Type_create_resized(MPI_HBTHalo_Id_t,(MPI_Aint)0, extent, &MPI_HBTHalo_Id_t);
MPI_Type_commit(&MPI_HBTHalo_Id_t);
#undef NumAttr
}
void DistributeHaloes(mpi::communicator &world, int root, vector <Halo_t> & InHalos, vector <Halo_t> & OutHalos, const ParticleSnapshot_t &snap)
/*distribute InHalos from root to around world. 
 *the destination of each halo is the one whose particle snapshot holds the most of this halo's particles. 
 *the distributed haloes are appended to OutHalos on each node.
*/
{
  typedef vector <Halo_t> HaloList_t;
  int thisrank=world.rank();
  HaloList_t HaloBuffer;
  HaloList_t & WorkingHalos=(world.rank()==root)?InHalos:HaloBuffer;
  MPI_Datatype MPI_HBT_Particle;
  create_MPI_Particle_type(MPI_HBT_Particle);
 
  //broadcast haloes
  {
  HBTInt nhalo=WorkingHalos.size();
  MPI_Bcast(&nhalo, 1, MPI_HBT_INT, root, world);
  if(world.rank()!=root)
	WorkingHalos.resize(nhalo);
  //send id
  MPI_Datatype MPI_HBTHalo_Id_t;
  create_MPI_Halo_Id_type(MPI_HBTHalo_Id_t);
  MPI_Bcast(WorkingHalos.data(), nhalo, MPI_HBTHalo_Id_t, root, world);
  MPI_Type_free(&MPI_HBTHalo_Id_t);
  
  //send sizes and prepare buffer
  vector <HBTInt> HaloSizes(nhalo);
  if(world.rank()==root)
  {
	for(HBTInt i=0;i<nhalo;i++)
	  HaloSizes[i]=WorkingHalos[i].Particles.size();
  }
  MPI_Bcast(HaloSizes.data(), nhalo, MPI_HBT_INT, root, world);
  if(world.rank()!=root)
  {
	for(HBTInt i=0;i<nhalo;i++)
	  WorkingHalos[i].Particles.resize(HaloSizes[i]);
  }
  
  //send particles
  vector <MPI_Aint> HaloBuffer(nhalo);
  for(HBTInt i=0;i<nhalo;i++)
  {
	 MPI_Aint p;
	 MPI_Address(WorkingHalos[i].Particles.data(),&p);
	 HaloBuffer[i]=p;
  }
  MPI_Datatype HaloType;
  MPI_Type_create_hindexed(nhalo, HaloSizes.data(), HaloBuffer.data(), MPI_HBT_Particle, &HaloType);
  MPI_Type_commit(&HaloType);
  MPI_Bcast(MPI_BOTTOM, 1, HaloType, root, world);
  MPI_Type_free(&HaloType);
  }
  
  //fill particles and decide movement
  struct SizeRank_t
  {
	HBTInt np;
	int rank;
  };
  vector <SizeRank_t> size(WorkingHalos.size()), maxsize(WorkingHalos.size());
  for(HBTInt haloid=0;haloid<WorkingHalos.size();haloid++)
  {
	Halo_t::ParticleList_t & Particles=WorkingHalos[haloid].Particles;
	size[haloid].np=0;
	size[haloid].rank=thisrank;
	HBTInt &np=size[haloid].np;
	for(HBTInt i=0;i<Particles.size();i++)
	{
	  HBTInt index=snap.GetIndex(Particles[i]);
	  if(index!=SpecialConst::NullParticleId)
		Particles[np++]=snap.Particles[index];
	}
	Particles.resize(np);
  }
#ifdef HBT_INT8
#define MPI_HBTPair MPI_LONG_INT
#else
#define MPI_HBTPair MPI_2INT
#endif
	MPI_Allreduce(size.data(), maxsize.data(), size.size(), MPI_HBTPair, MPI_MAXLOC, world);
	
	vector <vector <MPI_Aint> > SendBuffers(world.size()), ReceiveBuffers(world.size());
	vector <vector<int> > SendSizes(world.size()), ReceiveSizes(world.size());
	vector <MPI_Datatype> SendTypes(world.size()), ReceiveTypes(world.size());
  
	for(HBTInt haloid=0;haloid<WorkingHalos.size();haloid++)//packing
	{
	  int rank=maxsize[haloid].rank;
	  auto & Particles=WorkingHalos[haloid].Particles;
	  MPI_Aint p;
	  MPI_Address(Particles.data(),&p);
	  SendBuffers[rank].push_back(p);
	  SendSizes[rank].push_back(Particles.size());
	}
	for(int rank=0;rank<world.size();rank++)
	{
	  if(accumulate(SendSizes[rank].begin(), SendSizes[rank].end(), 0L)>INT_MAX)
		throw runtime_error("Error: sending more than INT_MAX particles around with MPI causes overflow. try increase the number of mpi threads.\n");
	  MPI_Type_create_hindexed(SendSizes[rank].size(), SendSizes[rank].data(), SendBuffers[rank].data(), MPI_HBT_Particle, &SendTypes[rank]);
	  MPI_Type_commit(&SendTypes[rank]);
	}

	VectorAllToAll(world, SendSizes, ReceiveSizes, MPI_INT);
	
	HBTInt NumNewHalos=ReceiveSizes[0].size();
	OutHalos.resize(OutHalos.size()+NumNewHalos);
	auto NewHalos=OutHalos.end()-NumNewHalos;
	for(int rank=0;rank<world.size();rank++)
	  ReceiveBuffers[rank].resize(NumNewHalos);
	for(HBTInt haloid=0;haloid<NumNewHalos;haloid++)
	{
	  HBTInt np=0;
	  for(int rank=0; rank<world.size(); rank++)
		np+=ReceiveSizes[rank][haloid];
	  auto &Particles=NewHalos[haloid].Particles;
	  Particles.resize(np);
	  np=0;
	  for(int rank=0; rank<world.size(); rank++)
	  {
		MPI_Aint p;
		MPI_Address(Particles.data()+np, &p);
		ReceiveBuffers[rank][haloid]=p;
		np+=ReceiveSizes[rank][haloid];
	  }
	}	
	for(int rank=0;rank<world.size();rank++)
	{
	  MPI_Type_create_hindexed(NumNewHalos, ReceiveSizes[rank].data(), ReceiveBuffers[rank].data(), MPI_HBT_Particle, &ReceiveTypes[rank]);
	  MPI_Type_commit(&ReceiveTypes[rank]);
	}
	
	vector <int> Counts(world.size(),1), Disps(world.size(),0);
	MPI_Alltoallw(MPI_BOTTOM, Counts.data(), Disps.data(), SendTypes.data(), MPI_BOTTOM, Counts.data(), Disps.data(), ReceiveTypes.data(), world);
	for(int rank=0;rank<world.size();rank++)
	{
	  MPI_Type_free(&SendTypes[rank]);
	  MPI_Type_free(&ReceiveTypes[rank]);
	}
	MPI_Type_free(&MPI_HBT_Particle);
	
	//copy other properties
	for(HBTInt haloid=0,i=0;haloid<WorkingHalos.size();haloid++)
	{
	  if(maxsize[haloid].rank==thisrank)
		NewHalos[i++].HaloId=WorkingHalos[haloid].HaloId;
	}
}
void HaloSnapshot_t::UpdateParticles(mpi::communicator &world, const ParticleSnapshot_t &snap)
{
  HaloList_t LocalHalos;
  for(int rank=0;rank<world.size();rank++)//one by one through the nodes
	DistributeHaloes(world, rank, Halos, LocalHalos, snap);
  Halos.swap(LocalHalos);
  
  TotNumberOfParticles=0;
  NumPartOfLargestHalo=0;
  for(auto &&h: LocalHalos)
  {
	HBTInt np=h.Particles.size();
	TotNumberOfParticles+=np;//local
	if(NumPartOfLargestHalo<np) NumPartOfLargestHalo=np;//local
  }
}
/* deprecated.
void HaloSnapshot_t::ParticleIdToIndex(const ParticleSnapshot_t& snapshot)
{
  chrono::high_resolution_clock::time_point time_begin, time_end;
#pragma omp master
  {
  ParticleSnapshot=&snapshot;
  }
#pragma omp for //maybe try collapse(2)? need to remove intermediate variables to enable this.
  for(HBTInt haloid=0;haloid<Halos.size();haloid++)
  {
	Halo_t::ParticleList_t & Particles=Halos[haloid].Particles;
	HBTInt np=Particles.size();
	for(HBTInt pid=0;pid<np;pid++)
	  Particles[pid]=snapshot.GetIndex(Particles[pid]);//this should be safe since its a const func
  }
}

void HaloSnapshot_t::ParticleIndexToId()
{
#pragma omp for
  for(HBTInt haloid=0;haloid<Halos.size();haloid++)
  {
	Halo_t::ParticleList_t &Particles=Halos[haloid].Particles;
	HBTInt nP=Halos[haloid].Particles.size();
	for(HBTInt pid=0;pid<nP;pid++)
	  Particles[pid]=ParticleSnapshot->GetId(Particles[pid]);
  }
#pragma omp single
  ParticleSnapshot=nullptr;
}
*/
void HaloSnapshot_t::AverageCoordinates()
{
#pragma omp for
  for(HBTInt i=0;i<Halos.size();i++)
  {
	AveragePosition(Halos[i].ComovingPosition, Halos[i].Particles.data(), Halos[i].Particles.size());
	AverageVelocity(Halos[i].PhysicalVelocity, Halos[i].Particles.data(), Halos[i].Particles.size());
  }
}

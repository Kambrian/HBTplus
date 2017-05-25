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

#include "snapshot.h"
#include "mymath.h"
#include "mpi_wrapper.h"

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
	const int root=0;
	if(newrank==root)
	{
	  MPI_Reduce(MPI_IN_PLACE, &IdMin, 1, MPI_HBT_INT, MPI_MIN, root, newcomm);
	  MPI_Reduce(MPI_IN_PLACE, &IdMax, 1, MPI_HBT_INT, MPI_MAX, root, newcomm);
	  flag_contig=(IdMax-IdMin+1==NumberOfParticlesOnAllNodes);
	  if(flag_contig) cout<<"Contiguous particle Ids."<<endl;
	}
	else
	{
	  MPI_Reduce(&IdMin, &IdMin, 1, MPI_HBT_INT, MPI_MIN, root, newcomm);
	  MPI_Reduce(&IdMax, &IdMax, 1, MPI_HBT_INT, MPI_MAX, root, newcomm);
	  IdMin=0;
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
  
  {
    HBTInt np=Particles.size();
    MPI_Allreduce(&np, &NumberOfParticlesOnAllNodes, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
  }
  
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
  
  sort(Particles.begin(), Particles.end(), CompParticleId);
  IdMin=Particles.front().Id; 
  IdMax=Particles.back().Id;
  ProcessIdRanges.resize(world.size()+1);
  MPI_Allgather(&IdMin, 1, MPI_HBT_INT, ProcessIdRanges.data(), 1, MPI_HBT_INT, world.Communicator);
  ProcessIdRanges.back()=IdMax+1;
  MPI_Bcast(&ProcessIdRanges.back(), 1, MPI_HBT_INT, world.size()-1, world.Communicator);
  
//   cout<<Particles.size()<<" particles received on node "<<world.rank()<<": IdRange=("<<IdMin<<","<<IdMax<<") "<<endl;
}
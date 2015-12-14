#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"
#include "gravity_tree.h"
/*
void MemberShipTable_t::SubIdToTrackId(const SubhaloList_t& Subhalos)
{
   for(HBTInt i=0;i<AllMembers.size();i++)
	AllMembers[i]=Subhalos[AllMembers[i]].TrackId;
}
void MemberShipTable_t::TrackIdToSubId(SubhaloList_t& Subhalos)
{
cout<<"Warning: TrackIdToSubId ToBe fully Implemented!\n";
// exit(1);
}
*/
  
void SubhaloSnapshot_t::BuildMPIDataType()
{
/*to create the struct data type for communication*/	
Subhalo_t p;
#define NumAttr 9
MPI_Datatype oldtypes[NumAttr];
int blockcounts[NumAttr];
MPI_Aint   offsets[NumAttr], origin,extent;

MPI_Get_address(&p,&origin);
MPI_Get_address((&p)+1,&extent);//to get the extent of s
extent-=origin;

int i=0;
#define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
RegisterAttr(TrackId, MPI_HBT_INT, 1)
RegisterAttr(Nbound, MPI_HBT_INT, 1)
RegisterAttr(HostHaloId, MPI_HBT_INT, 1)
RegisterAttr(Rank, MPI_HBT_INT, 1)
RegisterAttr(LastMaxMass, MPI_HBT_INT, 1)
RegisterAttr(SnapshotIndexOfLastMaxMass, MPI_HBT_INT, 1)
RegisterAttr(SnapshotIndexOfLastIsolation, MPI_HBT_INT, 1)
RegisterAttr(ComovingPosition[0], MPI_HBT_REAL, 3)
RegisterAttr(PhysicalVelocity[0], MPI_HBT_REAL, 3)
assert(offsets[i-1]-offsets[i-2]==sizeof(HBTReal)*3);//to make sure HBTxyz is stored locally.
#undef RegisterAttr
assert(i==NumAttr);

MPI_Type_create_struct(NumAttr,blockcounts,offsets,oldtypes, &MPI_HBT_SubhaloShell_t);
MPI_Type_create_resized(MPI_HBT_SubhaloShell_t,(MPI_Aint)0, extent, &MPI_HBT_SubhaloShell_t);
MPI_Type_commit(&MPI_HBT_SubhaloShell_t);
#undef NumAttr
}
void SubhaloSnapshot_t::UpdateParticles(mpi::communicator& world, const ParticleSnapshot_t& snapshot)
{
  SetEpoch(snapshot);
  SubhaloList_t LocalSubhalos;
  for(int rank=0;rank<world.size();rank++)//one by one through the nodes
	DistributeHaloes(world, rank, Subhalos, LocalSubhalos, snapshot, MPI_HBT_SubhaloShell_t);
  Subhalos.swap(LocalSubhalos);
}
/*
void SubhaloSnapshot_t::ParticleIdToIndex(const ParticleSnapshot_t& snapshot)
{//also bind to snapshot
#pragma omp single
  SnapshotPointer=&snapshot;
#pragma omp for
	for(HBTInt subid=0;subid<Subhalos.size();subid++)
	{
	  Subhalo_t::ParticleList_t & Particles=Subhalos[subid].Particles;
	  HBTInt nP=Particles.size();
	  for(HBTInt pid=0;pid<nP;pid++)
		Particles[pid]=snapshot.GetIndex(Particles[pid]);
	}
}
void SubhaloSnapshot_t::ParticleIndexToId()
{
#pragma omp for
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	Subhalo_t::ParticleList_t & Particles=Subhalos[subid].Particles;
	HBTInt nP=Particles.size();
	for(HBTInt pid=0;pid<nP;pid++)
	  Particles[pid]=SnapshotPointer->GetId(Particles[pid]);
  }
#pragma omp single
  SnapshotPointer=nullptr;
}
*/
void SubhaloSnapshot_t::AverageCoordinates()
{
#pragma omp for
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	int coresize=GetCoreSize(Subhalos[subid].Nbound);
	AveragePosition(Subhalos[subid].ComovingPosition, Subhalos[subid].Particles.data(), coresize);
	AverageVelocity(Subhalos[subid].PhysicalVelocity, Subhalos[subid].Particles.data(), Subhalos[subid].Nbound);
  }
}
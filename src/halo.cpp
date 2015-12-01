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

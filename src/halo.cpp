#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <glob.h>
#include <algorithm>
#include <chrono>

#include "mpi_wrapper.h"
#include "mymath.h"
#include "halo.h"
#include "particle_exchanger.h"

// #include <cstdio>
// #include <cstdlib>

void create_MPI_Halo_Id_type(MPI_Datatype &MPI_HBTHalo_Id_t)
{
/*to create the struct containing only haloid*/	
Halo_t p;
#define NumAttr 4
MPI_Datatype oldtypes[NumAttr];
int blockcounts[NumAttr];
MPI_Aint   offsets[NumAttr], origin,extent;

MPI_Get_address(&p,&origin);
MPI_Get_address((&p)+1,&extent);//to get the extent of s
extent-=origin;

int i=0;
#define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
RegisterAttr(HaloId, MPI_HBT_INT, 1)
RegisterAttr(ComovingAveragePosition[0], MPI_HBT_REAL, 3)
RegisterAttr(PhysicalAverageVelocity[0], MPI_HBT_REAL, 3)
RegisterAttr(Mass, MPI_HBT_REAL, 1)
// assert(offsets[i-1]-offsets[i-2]==sizeof(HBTReal)*3);//to make sure HBTxyz is stored locally.
#undef RegisterAttr
assert(i==NumAttr);

MPI_Type_create_struct(i,blockcounts,offsets,oldtypes, &MPI_HBTHalo_Id_t);//some padding is added automatically by MPI as well
MPI_Type_create_resized(MPI_HBTHalo_Id_t,(MPI_Aint)0, extent, &MPI_HBTHalo_Id_t);
MPI_Type_commit(&MPI_HBTHalo_Id_t);
#undef NumAttr
}

void Halo_t::AverageCoordinates()
{
  AveragePosition(ComovingAveragePosition, Particles.data(), Particles.size());
  Mass=AverageVelocity(PhysicalAverageVelocity, Particles.data(), Particles.size());
}

void HaloSnapshot_t::BuildMPIDataType()
{
  create_MPI_Halo_Id_type(MPI_HBT_HaloId_t);
}
void HaloSnapshot_t::UpdateParticles(MpiWorker_t &world, const ParticleSnapshot_t &snap)
{
  SetEpoch(snap);
  if(!HBTConfig.GroupLoadedFullParticle)
  {
  HaloList_t LocalHalos;
  snap.ExchangeHalos(world, Halos, LocalHalos, MPI_HBT_HaloId_t);
  Halos.swap(LocalHalos);
  }
  
  TotNumberOfParticles=0;
  NumPartOfLargestHalo=0;
  for(auto &&h: Halos)
  {
	HBTInt np=h.Particles.size();
	TotNumberOfParticles+=np;//local
	if(NumPartOfLargestHalo<np) NumPartOfLargestHalo=np;//local
  }
}

class HaloParticleKeyList_t: public KeyList_t <HBTInt, HBTInt>
{
  typedef HBTInt Index_t;
  typedef HBTInt Key_t;
  vector <HBTInt> ParticleIds;
  vector <HBTInt> HaloIds;//local haloid
public:
  HaloParticleKeyList_t(HaloSnapshot_t &snap)
  {
	ParticleIds.reserve(snap.TotNumberOfParticles);
	HaloIds.reserve(snap.TotNumberOfParticles);
	for(HBTInt i=0;i<snap.Halos.size();i++)
	{
	  auto &Part=snap.Halos[i].Particles;
	  for(auto && p: Part)
	  {
		ParticleIds.push_back(p.Id);
		HaloIds.push_back(i);//local haloid
	  }
	}
  };
  Index_t size() const
  {
	return ParticleIds.size();
  }
  Key_t GetKey(Index_t i) const
  {
	return ParticleIds[i];
  }
  Index_t GetIndex(Index_t i) const
  {
	return HaloIds[i];
  }
};

void HaloSnapshot_t::FillParticleHash()
{
  HaloParticleKeyList_t Ids(*this); 
  ParticleHash.Fill(Ids, SpecialConst::NullHaloId);
}
void HaloSnapshot_t::ClearParticleHash()
{
  ParticleHash.Clear();
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

#ifndef PARTICLE_EXCHANGER_H_INCLUDED
#define PARTICLE_EXCHANGER_H_INCLUDED

#include <numeric>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <list>
#include <forward_list>

#include "datatypes.h"
#include "mymath.h"
#include "mpi_wrapper.h"
#include "snapshot.h"
#include "halo_particle_iterator.h"
#include "subhalo.h"

class OrderedParticle_t: public Particle_t
{
public:
  HBTInt Order;
//   using Particle_t::Particle_t;
  using Particle_t::operator=;//inherit assignment operator
  OrderedParticle_t(){};
  OrderedParticle_t(HBTInt id, HBTInt order): Particle_t(id), Order(order)
  {
  }
  OrderedParticle_t(const Particle_t &p, HBTInt order): Particle_t(p), Order(order)
  {
  }
};

class RemoteParticle_t: public OrderedParticle_t
{
public:
  int ProcessorId;
//   using Particle_t::Particle_t;
  using Particle_t::operator=;//inherit assignment operator
  RemoteParticle_t(){};
  RemoteParticle_t(HBTInt id, HBTInt order): OrderedParticle_t(id, order)
  {
  }
  RemoteParticle_t(const Particle_t& p, HBTInt order): OrderedParticle_t(p, order)
  {
  }
//   using OrderedParticle_t::OrderedParticle_t; //not supported by old compilers
};

namespace ParticleExchangeComp
{
  inline bool CompParticleOrder(const OrderedParticle_t &a, const OrderedParticle_t &b)
  {
	return a.Order<b.Order;
  }
  inline bool CompParticleId(const Particle_t &a, const Particle_t &b)
  {
	return a.Id<b.Id;
  }

  inline bool CompIdAndOrder(const OrderedParticle_t &a, const OrderedParticle_t &b)
  {
	bool a_type=(a.Id!=SpecialConst::NullParticleId);
	bool b_type=(b.Id!=SpecialConst::NullParticleId);

	if(a_type>b_type)
	  return true;
	
	if(a_type&&b_type)
	  return a.Order<b.Order; 
  }
  
  template <class OrderedParticleList_T>
  void RestoreParticleOrder(OrderedParticleList_T &P)
  {
	for(HBTInt i=0;i<P.size();i++)
	{
	  auto &p=P[i];
	  auto &j=p.Order;
	  while(i!=j)
		swap(p, P[j]);
	}
  }
}

#include "hash_remote.tpp"

extern void create_Mpi_RemoteParticleType(MPI_Datatype& dtype);

template <class Halo_T>
class ParticleExchanger_t
{
  MPI_Datatype MPI_HBTParticle_t;
  MpiWorker_t &world;
  const ParticleSnapshot_t &snap;
  vector <Halo_T> &InHalos;
  
  vector <HBTInt> HaloSizes, HaloOffsets;
  
  vector <RemoteParticle_t> LocalParticles;
  vector <OrderedParticle_t> RoamParticles;
  typedef vector <RemoteParticle_t>::iterator LocalParticleIterator_t;
  typedef vector <OrderedParticle_t>::iterator RoamParticleIterator_t;
  vector <LocalParticleIterator_t> LocalIterators;
  vector <RoamParticleIterator_t> RoamIterators;
  vector <HBTInt> LocalSizes, RoamSizes;

  void CollectParticles();
  void SendParticles();
  void QueryParticles();
  void RecvParticles();
  void RestoreParticles();
public:
  ParticleExchanger_t(MpiWorker_t &world, const ParticleSnapshot_t &snap, vector <Halo_T> &InHalos);
  ~ParticleExchanger_t()
  {
	MPI_Type_free(&MPI_HBTParticle_t);
  }
  void Exchange();
};

template <class Halo_T>
ParticleExchanger_t<Halo_T>::ParticleExchanger_t(MpiWorker_t &world, const ParticleSnapshot_t &snap, vector <Halo_T> &InHalos): world(world),snap(snap), InHalos(InHalos), LocalParticles(), RoamParticles()
{
  Particle_t().create_MPI_type(MPI_HBTParticle_t);
}

template <class Halo_T>
void ParticleExchanger_t<Halo_T>::CollectParticles()
{
  HaloSizes.reserve(InHalos.size());
//   HaloOffsets.reserve(InHalos.size());
  HBTInt np=0;
  for(auto &&h: InHalos) 
  {
// 	HaloOffsets.push_back(np);
	HBTInt n=h.Particles.size();
	HaloSizes.push_back(n);
	np+=n;
  }
  LocalParticles.reserve(np);
  np=0;
  for(auto &&h: InHalos)
  {
	for(auto &&p: h.Particles)
	{
	  LocalParticles.emplace_back(move(p), np++);
	}
	h.Particles.clear();//or hard clear to save memory?
  }
}

template <class Halo_T>
void ParticleExchanger_t<Halo_T>::SendParticles()
{
  sort(LocalParticles.begin(), LocalParticles.end(), ParticleExchangeComp::CompParticleId);
  
  LocalSizes.resize(world.size());
  LocalIterators.resize(world.size());
  int rank=0;
  for(auto it=LocalParticles.begin(); it!=LocalParticles.end();++it)
  {
	while(it->Id>=snap.ProcessIdRanges[rank])
	{
	  if(rank==world.size()) //no particle id should exceed ProcessIdRange.back()
	  {    
	    cerr<<"Error: invalid particle id: "<<it->Id<<endl;
	    exit(1);
	  }
	  LocalIterators[rank++]=it;
	}
  }
  while(rank<world.size())
	LocalIterators[rank++]=LocalParticles.end();
  for(rank=0;rank<world.size()-1;rank++)
	LocalSizes[rank]=LocalIterators[rank+1]-LocalIterators[rank];
  LocalSizes[rank]=LocalParticles.end()-LocalIterators.back();
  
  RoamSizes.resize(world.size());
  MPI_Alltoall(LocalSizes.data(), 1, MPI_HBT_INT, RoamSizes.data(), 1, MPI_HBT_INT, world.Communicator);
  HBTInt np=accumulate(RoamSizes.begin(), RoamSizes.end(), (HBTInt)0);
  RoamParticles.resize(np);
  RoamIterators.resize(world.size());
  auto it=RoamParticles.begin();
  for(int i=0;i<world.size();i++)
  {
	RoamIterators[i]=it;
	it+=RoamSizes[i];
  }
  
  MyAllToAll<Particle_t, LocalParticleIterator_t, RoamParticleIterator_t>(world, LocalIterators, LocalSizes, RoamIterators, MPI_HBTParticle_t);
}

template <class Halo_T>
void ParticleExchanger_t<Halo_T>::QueryParticles()
{
  HBTInt order=0;
  for(auto &&p: RoamParticles)	p.Order=order++;
  sort(RoamParticles.begin(), RoamParticles.end(), ParticleExchangeComp::CompParticleId);
  
  snap.GetIndices(RoamParticles);
  for(auto &&p: RoamParticles)
  {
    if(p.Id!=SpecialConst::NullParticleId)
      p=snap.Particles[p.Id];
  }
  
  ParticleExchangeComp::RestoreParticleOrder(RoamParticles);
}

template <class Halo_T>
void ParticleExchanger_t<Halo_T>::RecvParticles()
{
  MyAllToAll<Particle_t, RoamParticleIterator_t, LocalParticleIterator_t>(world, RoamIterators, RoamSizes, LocalIterators, MPI_HBTParticle_t);
  vector <OrderedParticle_t>().swap(RoamParticles);
  
  for(int rank=0;rank<world.size();rank++)
  {
	auto it_end=LocalIterators[rank]+LocalSizes[rank];
	for(auto it=LocalIterators[rank];it!=it_end;++it)
	  it->ProcessorId=rank;
  }
  
  ParticleExchangeComp::RestoreParticleOrder(LocalParticles);
}

template <class Halo_T>
void ParticleExchanger_t<Halo_T>::RestoreParticles()
{
  auto it=LocalParticles.begin();
  for(HBTInt ihalo=0;ihalo<InHalos.size();ihalo++)
  {
	auto it_end=it+HaloSizes[ihalo];
	InHalos[ihalo].Particles.assign(it, it_end);
	it=it_end;
  }
  assert(it==LocalParticles.end());
  vector <RemoteParticle_t>().swap(LocalParticles);

  for(auto &&h: InHalos)
    h.KickNullParticles();
}

template <class Halo_T>
void ParticleExchanger_t<Halo_T>::Exchange()
{
  CollectParticles();
  SendParticles();
  QueryParticles();
  RecvParticles();
  RestoreParticles();
}

template <class Halo_T>
void DecideTargetProcessor(int NumProc, vector <Halo_T> &InHalos, vector <IdRank_t> &TargetRank)
{
  auto dims=ClosestFactors(NumProc, 3);
  HBTxyz step;
  for(int i=0;i<3;i++)
	step[i]=HBTConfig.BoxSize/dims[i];
  
#pragma omp parallel for
  for(HBTInt i=0;i<InHalos.size();i++)
  {
	InHalos[i].AverageCoordinates();
	TargetRank[i].Id=i;
	TargetRank[i].Rank=AssignCell(InHalos[i].ComovingAveragePosition, step, dims);
  }
}

template <class Halo_T>
void ParticleSnapshot_t::ExchangeHalos(MpiWorker_t& world, vector <Halo_T>& InHalos, vector<Halo_T>& OutHalos, MPI_Datatype MPI_Halo_Shell_Type) const
{
  typedef typename vector <Halo_T>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;
  
//   cout<<"Query particle..."<<flush;
  {//query particles
	ParticleExchanger_t <Halo_T>Exchanger(world, *this, InHalos);
	Exchanger.Exchange();
  }
  
  vector <IdRank_t>TargetRank(InHalos.size());
  DecideTargetProcessor(world.size(), InHalos, TargetRank);

  //distribute halo shells
	vector <int> SendHaloCounts(world.size(),0), RecvHaloCounts(world.size()), SendHaloDisps(world.size()), RecvHaloDisps(world.size());
	sort(TargetRank.begin(), TargetRank.end(), CompareRank);
	vector <Halo_T> InHalosSorted(InHalos.size());
	vector <HBTInt> InHaloSizes(InHalos.size());
	for(HBTInt haloid=0;haloid<InHalos.size();haloid++)
	{
	  InHalosSorted[haloid]=move(InHalos[TargetRank[haloid].Id]);
	  SendHaloCounts[TargetRank[haloid].Rank]++;
	  InHaloSizes[haloid]=InHalosSorted[haloid].Particles.size();
	}
	MPI_Alltoall(SendHaloCounts.data(), 1, MPI_INT, RecvHaloCounts.data(), 1, MPI_INT, world.Communicator);
	CompileOffsets(SendHaloCounts, SendHaloDisps);
	HBTInt NumNewHalos=CompileOffsets(RecvHaloCounts, RecvHaloDisps);
	OutHalos.resize(OutHalos.size()+NumNewHalos);
	auto NewHalos=OutHalos.end()-NumNewHalos;
	MPI_Alltoallv(InHalosSorted.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_Halo_Shell_Type, &NewHalos[0], RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_Halo_Shell_Type, world.Communicator);
  //resize receivehalos
	vector <HBTInt> OutHaloSizes(NumNewHalos);
	MPI_Alltoallv(InHaloSizes.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HBT_INT, OutHaloSizes.data(), RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HBT_INT, world.Communicator);
	for(HBTInt i=0;i<NumNewHalos;i++)
	  NewHalos[i].Particles.resize(OutHaloSizes[i]);
	
	{
	//distribute halo particles
	MPI_Datatype MPI_HBT_Particle;
	Particle_t().create_MPI_type(MPI_HBT_Particle);
	//create combined iterator for each bunch of haloes
	vector <ParticleIterator_t> InParticleIterator(world.size());
	vector <ParticleIterator_t> OutParticleIterator(world.size());
	for(int rank=0;rank<world.size();rank++)
	{
	  InParticleIterator[rank].init(InHalosSorted.begin()+SendHaloDisps[rank], InHalosSorted.begin()+SendHaloDisps[rank]+SendHaloCounts[rank]);
	  OutParticleIterator[rank].init(NewHalos+RecvHaloDisps[rank], NewHalos+RecvHaloDisps[rank]+RecvHaloCounts[rank]);
	}
	vector <HBTInt> InParticleCount(world.size(),0);
	for(HBTInt i=0;i<InHalosSorted.size();i++)
	  InParticleCount[TargetRank[i].Rank]+=InHalosSorted[i].Particles.size();
	
	MyAllToAll<Particle_t, ParticleIterator_t, ParticleIterator_t>(world, InParticleIterator, InParticleCount, OutParticleIterator, MPI_HBT_Particle);
	
	MPI_Type_free(&MPI_HBT_Particle);
	}
}
#endif
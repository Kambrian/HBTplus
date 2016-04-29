#ifndef PARTICLE_EXCHANGER_H_INCLUDED
#define PARTICLE_EXCHANGER_H_INCLUDED

#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <list>
#include <forward_list>

#include "datatypes.h"
#include "mymath.h"
#include "mpi_wrapper.h"

struct RemoteParticleId_t
{
  HBTInt Id;
  int ProcessorId;
  HBTInt Order;
  RemoteParticleId_t()=default;
  RemoteParticleId_t(HBTInt id, int processorId, HBTInt order): Id(id), ProcessorId(processorId), Order(order)
  {
  }
};
class RemoteParticle_t: public Particle_t
{
public:
  int ProcessorId;
  HBTInt Order;
  
//   using Particle_t::Particle_t;
  using Particle_t::operator=;//inherit assignment operator
  RemoteParticle_t()=default;
  RemoteParticle_t(HBTInt id, int processorId, HBTInt order): Particle_t(id), ProcessorId(processorId), Order(order)
  {
  }
  RemoteParticle_t(const Particle_t &p, int processorId, HBTInt order): Particle_t(p), ProcessorId(processorId), Order(order)
  {
  }
};

template <class Pair_t, class Val_t>
inline int CompPairWithValue(const Pair_t a, const Val_t b)
{
  return (a.Key<b);
};
template <class Key_t, class Index_t>
void MappedIndexTable_t<Key_t, Index_t>::GetIndices(ParticleIdList_T &particles) const
{ 
#define ALWAYS_BATCH_BINARY_SEARCH
  
#ifdef ALWAYS_BATCH_BINARY_SEARCH
  GetIndicesRecursive(particles, 0, particles.size(), Map.begin(), Map.end());//batch-binary-search: is this always faster?
#else  
  if(particles.size()<NumQueryCrit)//do individual binary search
  {
	for(auto &&p: particles)
	  p.Id=GetIndex(p.Id);
	return;
  }
  //otherwise do batch sequential search
  auto &null=BaseClass_t::NullIndex;
  
  auto it_p=particles.begin();
  auto it_map=Map.begin();
  while(true)
  {
	if(it_p==particles.end()) return;
	
	if(it_map==Map.end()) break;
	
	if(it_p->Id<it_map->Key)
	{
	  it_p->Id=null;
	  ++it_p;
	}
	else if(it_p->Id==it_map->Key)
	{
	  it_p->Id=it_map->Index;
	  ++it_p;
	}
	else
	  ++it_map;
  }
  
  while(true)
  {
	  it_p->Id=null;
	  ++it_p;
	  if(it_p==particles.end()) return;
  }
#endif  
}

template <class Key_t, class Index_t>
void MappedIndexTable_t<Key_t, Index_t>::GetIndicesRecursive(ParticleIdList_T &particles, HBTInt imin, HBTInt imax, MapIter_t MapBegin,  MapIter_t MapEnd) const
{
  //GetIndices of particles in storage range [imin, imax) from map [MapBegin, MapEnd).
  auto &null=BaseClass_t::NullIndex;
 
  if(MapBegin==MapEnd)
  {
	for(HBTInt i=imin;i<imax;i++)
	  particles[i].Id=null;
	return;
  }
  
  if(imin>=imax) return;
  
  HBTInt imid;
  if(imax-imin==1) 
	imid=imin;
  else
	imid=(imin+imax)/2;
  Key_t key=particles[imid].Id;
  MapIter_t MapMid=lower_bound(MapBegin, MapEnd, key, CompPairWithValue<Pair_t, Key_t>);
  MapIter_t MapEndLeft=MapMid, MapBeginRight=MapMid;
  if(MapMid==MapEnd||MapMid->Key>key)
	particles[imid].Id=null;
  else
  {
	particles[imid].Id=MapMid->Index;
	++MapEndLeft;
  }
	
  GetIndicesRecursive(particles, imin, imid, MapBegin, MapEndLeft);
  GetIndicesRecursive(particles, imid+1, imax, MapBeginRight, MapEnd);
}

template <class Key_t, class Index_t>
void FlatIndexTable_t<Key_t, Index_t>::GetIndices(ParticleIdList_T &particles) const
{
  for(auto &&p: particles)
	p.Id=GetIndex(p.Id);
}


extern void create_Mpi_RemoteParticleType(MPI_Datatype& dtype);
extern void create_Mpi_RemoteParticleIdType(MPI_Datatype& dtype);
class ParticleExchanger_t
{
  MPI_Datatype MPI_RemoteParticleId_t, MPI_RemoteParticle_t;//todo: init and free them
  MpiWorker_t &world;
  const ParticleSnapshot_t &snap;
  typedef vector <RemoteParticleId_t> ParticleStack_t;
  ParticleStack_t ParticlesToProcess, ParticlesToSend;
  HBTInt SendStackSize, SendStackSize0;
  int CurrSendingRank;
  vector <RemoteParticle_t> LocalParticles;
public:
  template <class HaloParticleIterator_t>
  ParticleExchanger_t(MpiWorker_t &_world, const ParticleSnapshot_t &_snap, HaloParticleIterator_t &particle_it);
  ~ParticleExchanger_t()
  {
	MPI_Type_free(&MPI_RemoteParticleId_t);
	MPI_Type_free(&MPI_RemoteParticle_t);
  }
  void BcastParticles(HBTInt &ParticleCount);
  bool GatherParticles(HBTInt capacity);
  void QueryParticles();
  void Exchange();
  void RestoreParticles();
  template <class Halo_T>
  void UnPackHaloParticles(vector <Halo_T> &InHalos, vector <IdRank_t> &TargetRank);
};

template <class HaloParticleIterator_t>
ParticleExchanger_t::ParticleExchanger_t(MpiWorker_t &_world, const ParticleSnapshot_t &_snap, HaloParticleIterator_t &particle_it): world(_world), snap(_snap),  CurrSendingRank(0)
{
  create_Mpi_RemoteParticleType(MPI_RemoteParticle_t);
  create_Mpi_RemoteParticleIdType(MPI_RemoteParticleId_t);
  {
	SendStackSize=0;
	while(!particle_it.is_end())
	{
	  ++particle_it;
	  ++SendStackSize;
	}
	particle_it.reset();
	ParticlesToSend.reserve(SendStackSize);
  }
  HBTInt order=0;
  while(!particle_it.is_end())
  {
	ParticlesToSend.emplace_back((*particle_it).Id, world.rank(), order);
	++particle_it;
	++order;
  }
  SendStackSize=order;
  SendStackSize0=SendStackSize;//backup size for later sanity check
  particle_it.reset();
}

template <class Halo_T>
void ParticleExchanger_t::UnPackHaloParticles(vector <Halo_T> &InHalos, vector <IdRank_t> &TargetRank)
{
  TargetRank.reserve(InHalos.size());
  auto it=LocalParticles.begin();
  for(HBTInt ihalo=0;ihalo<InHalos.size();ihalo++)
  {
	auto &h=InHalos[ihalo];
	vector <HBTInt> counter(world.size(),0);
	for(auto &&p: h.Particles)
	{
	  auto &rp=*it;
	  p=rp;
	  if(rp.ProcessorId>=0)
		counter[rp.ProcessorId]++;
	  ++it;
	}
	int targetrank=world.rank();
	auto target=max_element(counter.begin(), counter.end());
	if(*target) targetrank=target-counter.begin();
	TargetRank.emplace_back(ihalo, targetrank);
  }
  vector <RemoteParticle_t>().swap(LocalParticles);//clear up
}

template <class HaloIterator>
class HaloParticleIterator_t
{
  typedef vector<Particle_t>::iterator particle_iterator;
  HaloIterator FirstHalo, EndHalo, CurrHalo;
  particle_iterator CurrPart;
public:
  HaloParticleIterator_t()=default;
  HaloParticleIterator_t(const HaloIterator &begin, const HaloIterator &end)
  {
	init(begin, end);
  }
  void init(HaloIterator begin, HaloIterator end)
  {
	while((begin!=end)&&(begin->Particles.size()==0))//skip empty ones, though not necessary for current HBT2
	  ++begin;
	FirstHalo=begin;
	EndHalo=end;
	reset();
  }
  void reset()
  {
	CurrHalo=FirstHalo;
	if(CurrHalo!=EndHalo)
	  CurrPart=FirstHalo->Particles.begin();
  }
  particle_iterator begin()
  {
	return FirstHalo->Particles.begin();
  }
  HaloParticleIterator_t<HaloIterator> & operator ++()//left operator
  {
	++CurrPart;
	while(CurrPart==CurrHalo->Particles.end())//increment halo and skip empty haloes
	{
	  ++CurrHalo;
	  if(CurrHalo==EndHalo) break;
	  CurrPart=CurrHalo->Particles.begin();
	}
	return *this;
  }
  Particle_t & operator *()
  {
	return *CurrPart;
  }
  bool is_end()
  {
	return CurrHalo==EndHalo;
  }
};

template <class ParticleIdList_t>
void ParticleSnapshot_t::GetIndices(ParticleIdList_t& particles) const
{//ParticleIdList_t is a list of particle structs containing at least an Id field
  return ParticleHash->GetIndices(particles);
}

template <class Halo_T>
void ParticleSnapshot_t::ExchangeHalos(MpiWorker_t& world, vector <Halo_T>& InHalos, vector<Halo_T>& OutHalos, MPI_Datatype MPI_Halo_Shell_Type) const
{
  typedef typename vector <Halo_T>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;
  
//   cout<<"Query particle..."<<flush;
  vector <IdRank_t>TargetRank;
  {//query particles
	ParticleIterator_t InParticles(InHalos.begin(), InHalos.end());
	ParticleExchanger_t Exchanger(world, *this, InParticles);
	Exchanger.Exchange();
	Exchanger.UnPackHaloParticles(InHalos, TargetRank);
  }
  
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
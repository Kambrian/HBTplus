#ifndef PARTICLE_EXCHANGER_H_INCLUDED
#define PARTICLE_EXCHANGER_H_INCLUDED

#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include <list>
#include <forward_list>

#include "datatypes.h"
#include "mymath.h"
#include "mpi_wrapper.h"

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
  RemoteParticle_t(Particle_t &p, int processorId, HBTInt order): Particle_t(p), ProcessorId(processorId), Order(order)
  {
  }
};

extern void create_Mpi_RemoteParticleType(MPI_Datatype& dtype, bool IdOnly=false);
class ParticleExchanger_t
{
  int PrevRank, NextRank;
  int iloop, nloop;
  const HBTInt EndParticleId;
  const int TagQuery;
  const int maxbuffersize;
  vector <RemoteParticle_t> SendBuffer, RecvBuffer;
  MPI_Request ReqSend;
  int ChannelIsClean;
  MPI_Datatype MPI_RemoteParticleId_t, MPI_RemoteParticle_t;//todo: init and free them
  MpiWorker_t &world;
  const ParticleSnapshot_t &snap;
  typedef list <RemoteParticle_t> ParticleStack_t;
  ParticleStack_t ParticlesToProcess;
  ParticleStack_t ParticlesToSend;
  vector <RemoteParticle_t> LocalParticles;
public:
  template <class HaloParticleIterator_t>
  ParticleExchanger_t(MpiWorker_t &_world, const ParticleSnapshot_t &_snap, HaloParticleIterator_t &particle_it): world(_world), snap(_snap), iloop(0), EndParticleId(-1), TagQuery(1), maxbuffersize(100), ReqSend(MPI_REQUEST_NULL), ChannelIsClean(1), SendBuffer(maxbuffersize), RecvBuffer(maxbuffersize)
  {
	nloop=world.size();
	PrevRank=world.prev();
	NextRank=world.next();
	create_Mpi_RemoteParticleType(MPI_RemoteParticle_t);
	create_Mpi_RemoteParticleType(MPI_RemoteParticleId_t, true);
	HBTInt order=0;
	while(!particle_it.is_end())
	{
	  ParticlesToProcess.emplace_back(*particle_it, world.rank(), order);
	  ++particle_it;
	  ++order;
	}
	particle_it.reset();
	ParticlesToProcess.emplace_back(EndParticleId, world.rank(), -1);
  }
  ~ParticleExchanger_t()
  {
	MPI_Type_free(&MPI_RemoteParticleId_t);
	MPI_Type_free(&MPI_RemoteParticle_t);
  }
  bool ProcessParticle(ParticleStack_t::iterator &it);
  void Exchange();
  void SendParticles();
  void ReceiveParticles(bool blocking);
  void RestoreParticles();
  bool FlushChannel()
  {
	if(!ChannelIsClean)
	  MPI_Test(&ReqSend, &ChannelIsClean, MPI_STATUS_IGNORE);
  }
  void WaitChannel()
  {
	if(!ChannelIsClean)
	  MPI_Wait(&ReqSend, MPI_STATUS_IGNORE);
  }
  bool AllDone()
  {
	return iloop>=nloop;
  }
  template <class Halo_T>
  void UnPackHaloParticles(vector <Halo_T> &InHalos, vector <IdRank_t> &TargetRank);
};
template <class Halo_T>
void ParticleExchanger_t::UnPackHaloParticles(vector <Halo_T> &InHalos, vector <IdRank_t> &TargetRank)
{
  TargetRank.reserve(InHalos.size());
  auto it=LocalParticles.begin();
  for(HBTInt ihalo=0;ihalo<InHalos.size();ihalo++)
  {
	auto &h=InHalos[ihalo];
	unordered_map<int, HBTInt> counter;
	for(auto &&p: h.Particles)
	{
	  auto &rp=*it;
	  p=rp;
	  if(rp.ProcessorId>=0)
		counter[rp.ProcessorId]++;
	  ++it;
	}
	int targetrank=world.rank(), n=0;
	for(auto &&x: counter)
	{
	  if(x.second>n)
	  {
		n=x.second;
		targetrank=x.first;
	  }
	}
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
  
//   cout<<"sending shells...\n";
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
	
// 	cout<<"sending particles...";
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
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

class RemoteParticle_t: public Particle_t
{
public:
  HBTInt Order;
  
//   using Particle_t::Particle_t;
  using Particle_t::operator=;//inherit assignment operator
  RemoteParticle_t()=default;
  RemoteParticle_t(HBTInt id, HBTInt order): Particle_t(id), Order(order)
  {
  }
  RemoteParticle_t(const Particle_t &p, HBTInt order): Particle_t(p), Order(order)
  {
  }
};

namespace ParticleExchangeComp
{
  inline bool CompParticleOrder(const RemoteParticle_t &a, const RemoteParticle_t &b)
  {
	return a.Order<b.Order;
  }
  inline bool CompParticleId(const RemoteParticle_t &a, const RemoteParticle_t &b)
  {
	return a.Id<b.Id;
  }

  inline bool CompIdAndOrder(const RemoteParticle_t &a, const RemoteParticle_t &b)
  {
	bool a_type=(a.Id!=SpecialConst::NullParticleId);
	bool b_type=(b.Id!=SpecialConst::NullParticleId);

	if(a_type>b_type)
	  return true;
	
	if(a_type&&b_type)
	  return a.Order<b.Order; 
  }
  
  extern void SortRemoteParticles(vector <RemoteParticle_t> &P);
  extern void ReduceRemoteParticle( void *in, void *inout, int *len, MPI_Datatype *dptr );
}

#include "hash_remote.tpp"

extern void create_Mpi_RemoteParticleType(MPI_Datatype& dtype);

template <class HaloParticleIterator_T>
class ParticleExchanger_t
{
  MPI_Op MPI_ReduceRemoteParticleOp;
  MPI_Datatype MPI_RemoteParticle_t;//todo: init and free them
  MpiWorker_t &world;
  const ParticleSnapshot_t &snap;
  vector <RemoteParticle_t> ParticlesToProcess;
  HBTInt NumPartToSend, NumPartToRecv;
  int CurrSendingRank, CurrRecvingRank;
  HBTInt CurrSendingOrder;
  HaloParticleIterator_T SendParticleIter, RecvParticleIter;
  vector <int> TargetProcessor;
  vector <int>::iterator TargetProcessorIter;
public:
  ParticleExchanger_t(MpiWorker_t &_world, const ParticleSnapshot_t &_snap, HaloParticleIterator_T &particle_it);
  ~ParticleExchanger_t()
  {
	MPI_Type_free(&MPI_RemoteParticle_t);
	MPI_Op_free(&MPI_ReduceRemoteParticleOp);
  }
  void BcastParticles(HBTInt &ParticleCount);
  bool GatherParticles(HBTInt capacity);
  void QueryParticles();
  void ReduceParticles(HBTInt &ParticleOffset, HBTInt &ParticleCount);
  bool RestoreParticles();
  void Exchange();
  template <class Halo_T>
  void CompileTargets(vector <Halo_T> &InHalos, vector <IdRank_t> &TargetRank);
};

template <class HaloParticleIterator_T>
ParticleExchanger_t<HaloParticleIterator_T>::ParticleExchanger_t(MpiWorker_t &_world, const ParticleSnapshot_t &_snap, HaloParticleIterator_T &particle_it): SendParticleIter(particle_it), RecvParticleIter(particle_it), world(_world), snap(_snap),  CurrSendingRank(0), CurrSendingOrder(0), CurrRecvingRank(0)
{
  create_Mpi_RemoteParticleType(MPI_RemoteParticle_t);
  MPI_Op_create( ParticleExchangeComp::ReduceRemoteParticle, true, &MPI_ReduceRemoteParticleOp ); 
  {
	NumPartToSend=0;
	while(!SendParticleIter.is_end())
	{
	  ++SendParticleIter;
	  ++NumPartToSend;
	}
	SendParticleIter.reset();
  }
  NumPartToRecv=NumPartToSend;//backup size for later sanity check
  TargetProcessor.resize(NumPartToSend);
  TargetProcessorIter=TargetProcessor.begin();
}

template <class HaloParticleIterator_T>
template <class Halo_T>
void ParticleExchanger_t<HaloParticleIterator_T>::CompileTargets(vector <Halo_T> &InHalos, vector <IdRank_t> &TargetRank)
{
  TargetRank.reserve(InHalos.size());
  auto it=TargetProcessor.begin();
  for(HBTInt ihalo=0;ihalo<InHalos.size();ihalo++)
  {
	auto &h=InHalos[ihalo];
	vector <HBTInt> counter(world.size(),0);
	for(auto &&p: h.Particles)
	{
	  if(*it>=0)
		counter[*it]++;
	  ++it;
	}
	int targetrank=world.rank();
	auto target=max_element(counter.begin(), counter.end());
	if(*target) targetrank=target-counter.begin();
	TargetRank.emplace_back(ihalo, targetrank);
  }
  vector <int>().swap(TargetProcessor);//clear up
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
	ParticleExchanger_t<ParticleIterator_t> Exchanger(world, *this, InParticles);
	Exchanger.Exchange();
	Exchanger.CompileTargets(InHalos, TargetRank);
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

template <class HaloParticleIterator_T>
void ParticleExchanger_t<HaloParticleIterator_T>::BcastParticles(HBTInt& ParticleCount)
{
  int root=CurrSendingRank;
  HBTInt np=ParticlesToProcess.size();
  MPI_Bcast(&ParticleCount, 1, MPI_HBT_INT, root, world.Communicator);
  if(ParticleCount)
  {
	//determine loops
	const int chunksize=1024*1024;
	HBTInt  Nloop=ceil(1.*ParticleCount/chunksize);
	int buffersize=ParticleCount/Nloop+1, nremainder=ParticleCount%Nloop;
	//transmit
	vector <HBTInt> buffer(buffersize);
	for(HBTInt iloop=0;iloop<Nloop;iloop++)
	{
	  if(iloop==nremainder)//switch sendcount from n+1 to n
	  {
		buffersize--;
		buffer.resize(buffersize);
	  }
	  if(world.rank()==root)//pack
	  {
		for(auto it_buff=buffer.begin();it_buff!=buffer.end();++it_buff)
		{
		  *it_buff=(*SendParticleIter).Id;
		  ++SendParticleIter;
		}
	  }
	  MPI_Bcast(buffer.data(), buffersize, MPI_HBT_INT, root, world.Communicator);
	  for(auto it_buff=buffer.begin();it_buff!=buffer.end();++it_buff)//unpack
	  {
		ParticlesToProcess.emplace_back(*it_buff, np++);
	  }
	}
  }
  if(world.rank()==root)
  {
	NumPartToSend-=ParticleCount;
	if(SendParticleIter.is_end()) CurrSendingRank++;
  }
  MPI_Bcast(&CurrSendingRank, 1, MPI_HBT_INT, root, world.Communicator);
}

template <class HaloParticleIterator_T>
bool ParticleExchanger_t<HaloParticleIterator_T>::GatherParticles(HBTInt capacity)
{
  HBTInt nsend;
  while(capacity)
  {
	if(world.rank()==CurrSendingRank)
	  nsend=min(NumPartToSend, capacity);
	BcastParticles(nsend);
	capacity-=nsend;
	if(CurrSendingRank==world.size()) return true;
  }
  
  return false;
}

template <class HaloParticleIterator_T>
void ParticleExchanger_t<HaloParticleIterator_T>::QueryParticles()
{
  sort(ParticlesToProcess.begin(), ParticlesToProcess.end(), ParticleExchangeComp::CompParticleId);

  snap.GetIndices(ParticlesToProcess);
  for(auto &&p: ParticlesToProcess)
	if(p.Id!=SpecialConst::NullParticleId)
	  p=snap.Particles[p.Id];
  
  ParticleExchangeComp::SortRemoteParticles(ParticlesToProcess);
//   sort(ParticlesToProcess.begin(), ParticlesToProcess.end(), ParticleExchangeComp::CompParticleOrder);
  
  HBTInt rank=world.rank();
  for(auto &&p: ParticlesToProcess)	p.Order=rank;
  
//  ParticlesToProcess.clear();
}

template <class HaloParticleIterator_T>
void ParticleExchanger_t<HaloParticleIterator_T>::Exchange()
{
  HBTInt capacity=ceil(1.*snap.NumberOfParticlesOnAllNodes/world.size());//decrease this if out of memory; increase this to increase efficiency
  while(true)
  {
	ParticlesToProcess.reserve(capacity);
	bool flag_end=GatherParticles(capacity);
	assert(ParticlesToProcess.size()<=capacity);
	QueryParticles();
	RestoreParticles();
	ParticlesToProcess.clear();
	if(flag_end) break;
  }
}

template <class HaloParticleIterator_T>
void ParticleExchanger_t<HaloParticleIterator_T>::ReduceParticles(HBTInt &ParticleOffset, HBTInt &ParticleCount)
{
  int root=CurrRecvingRank;
  MPI_Bcast(&ParticleCount, 1, MPI_HBT_INT, root, world.Communicator);
  RemoteParticle_t * data=ParticlesToProcess.data()+ParticleOffset;
  if(ParticleCount)
  {
	//determine loops
	const int chunksize=1024*1024;
	HBTInt  Nloop=ceil(1.*ParticleCount/chunksize);
	int buffersize=ParticleCount/Nloop+1, nremainder=ParticleCount%Nloop;
	//transmit
	HBTInt ndone=0;
	for(HBTInt iloop=0;iloop<Nloop;iloop++)
	{
	  if(iloop==nremainder)//switch sendcount from n+1 to n
		buffersize--;
	  if(world.rank()==root)
		MPI_Reduce(MPI_IN_PLACE, data+ndone, buffersize, MPI_RemoteParticle_t, MPI_ReduceRemoteParticleOp, root, world.Communicator);
	  else
		MPI_Reduce(data+ndone, NULL, buffersize, MPI_RemoteParticle_t, MPI_ReduceRemoteParticleOp, root, world.Communicator);
	  ndone+=buffersize;
	}
	if(world.rank()==root)//unload
	{
	  for(auto it_buff=data;it_buff!=data+ParticleCount;++it_buff)
	  {
		*TargetProcessorIter=it_buff->Order;
		++TargetProcessorIter;
		*RecvParticleIter=move(*it_buff);
		++RecvParticleIter;
	  }
	}
	ParticleOffset+=ParticleCount;
  }
  if(world.rank()==root)
  {
	NumPartToRecv-=ParticleCount;
	if(RecvParticleIter.is_end()) CurrRecvingRank++;
  }
  MPI_Bcast(&CurrRecvingRank, 1, MPI_HBT_INT, root, world.Communicator);
}

template <class HaloParticleIterator_T>
bool ParticleExchanger_t<HaloParticleIterator_T>::RestoreParticles()
{
  HBTInt nsend, offset=0;
  HBTInt capacity=ParticlesToProcess.size();
  while(capacity)
  {
	if(world.rank()==CurrRecvingRank)
	  nsend=min(NumPartToRecv, capacity);
	ReduceParticles(offset, nsend);
	capacity-=nsend;
	if(CurrRecvingRank==world.size()) return true;
  }
  
  return false;
}

#endif
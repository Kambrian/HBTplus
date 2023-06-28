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

struct HaloSlicer_t
{
    HBTInt ihalo_begin, ihalo_back, ipart_begin, ipart_end, np;//range of exchanging halos and particles in halos
    int BufferSize;//max number of particles to slice
    HaloSlicer_t(int buffersize):ihalo_begin(0), ihalo_back(0), ipart_begin(0), ipart_end(0), np(0), BufferSize(buffersize)
    {
      if(BufferSize==0) BufferSize=HBTConfig.ParticleExchangerBufferSize;
    }
    template <class Halo_T>
    int NextBuffer(const vector <Halo_T> &InHalos)
    {//select halo range to fill into buffer. return tot number of selected particles.
        np=0;
        ihalo_begin=ihalo_back;
        ipart_begin=ipart_end;
        if(ihalo_begin>=InHalos.size())
          return 0;
        if(ipart_begin>=InHalos[ihalo_begin].Particles.size())//overflow
        {
          ihalo_begin++;
          if(ihalo_begin>=InHalos.size())//end
            return 0;
          ipart_begin=0;
        }
        HBTInt ihalo=ihalo_begin;
        while(true)
        {
            HBTInt n_this=InHalos[ihalo].Particles.size();
            if(ihalo==ihalo_begin) n_this-=ipart_begin;
            np+=n_this;
            if(np>=BufferSize)//buffer full
            {
                HBTInt np_excess=np-BufferSize;
                ihalo_back=ihalo;
                ipart_end=InHalos[ihalo].Particles.size()-np_excess;
                np=BufferSize;
                break;
            }
            ihalo++;
            if(ihalo>=InHalos.size())//exceed
            {
              ihalo_back=InHalos.size()-1;
              ipart_end=InHalos.back().Particles.size();
              break;
            }
        }
        return np;
    }
    template <class Halo_T>
    void GetPartRange(const vector <Halo_T> &InHalos, HBTInt ihalo, HBTInt &i_begin, HBTInt &i_end)
    {//get particle range for ihalo
      i_begin=0;
      if(ihalo==ihalo_begin)
        i_begin=ipart_begin;
      i_end=InHalos[ihalo].Particles.size();
      if(ihalo==ihalo_back)
        i_end=ipart_end;
    }
};

template <class Halo_T>
class ParticleExchanger_t
{
  MPI_Datatype MPI_HBTParticle_t;
  MpiWorker_t &world;
  const ParticleSnapshot_t &snap;
  vector <Halo_T> &InHalos;
  HaloSlicer_t Slicer;
  vector <HBTInt> HaloSizes;


  vector <RemoteParticle_t> LocalParticles;
  vector <OrderedParticle_t> RoamParticles;
  typedef vector <RemoteParticle_t>::iterator LocalParticleIterator_t;
  typedef vector <OrderedParticle_t>::iterator RoamParticleIterator_t;
  vector <LocalParticleIterator_t> LocalIterators;
  vector <RoamParticleIterator_t> RoamIterators;
  vector <HBTInt> LocalSizes, RoamSizes;

  int CollectParticles();
  void SendParticles();
  void QueryParticles();
  void RecvParticles();
  void RestoreParticles();
public:
  ParticleExchanger_t(MpiWorker_t &world, const ParticleSnapshot_t &snap, vector <Halo_T> &InHalos, int buffer_size=0);//buffer_size=0 means using HBTConfig.ParticleExchangerBufferSize
  ~ParticleExchanger_t()
  {
	My_Type_free(&MPI_HBTParticle_t);
  }
  void Exchange();
};

template <class Halo_T>
ParticleExchanger_t<Halo_T>::ParticleExchanger_t(MpiWorker_t &world, const ParticleSnapshot_t &snap, vector <Halo_T> &InHalos, int buffer_size): world(world),snap(snap), InHalos(InHalos), Slicer(buffer_size), HaloSizes(),  LocalParticles(), RoamParticles()
{
  Particle_t().create_MPI_type(MPI_HBTParticle_t);
}

template <class Halo_T>
int ParticleExchanger_t<Halo_T>::CollectParticles()
{
  int np=Slicer.NextBuffer(InHalos);
  if(0==np)
    return 0;
  LocalParticles.clear();//not necessary since LocalParticles has been cleared during RestoreParticles.
  LocalParticles.reserve(np);
  int n=0;
  for(HBTInt ihalo=Slicer.ihalo_begin;ihalo<=Slicer.ihalo_back;ihalo++)
  {
    auto &h=InHalos[ihalo];
    HBTInt ipart, ipart_end;
    Slicer.GetPartRange(InHalos, ihalo, ipart, ipart_end);
    for(;ipart<ipart_end;ipart++)
      LocalParticles.emplace_back(move(h.Particles[ipart]), n++);
//     h.Particles.clear();//do not clear since we will need to know the particle sizes later
  }
  assert(n==np);
  return np;
}

template <class Halo_T>
void ParticleExchanger_t<Halo_T>::SendParticles()
{//order the particles and send them to the corresponding processing according to ProcessIdRange
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
      p=snap.Particles[p.Id];//query the particle property; may need to spawn particles due to star formation here.
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
  for(HBTInt ihalo=Slicer.ihalo_begin;ihalo<=Slicer.ihalo_back;ihalo++)
  {
    HBTInt ipart, ipart_end;
    Slicer.GetPartRange(InHalos, ihalo, ipart, ipart_end);
    for(;ipart<ipart_end;ipart++)
      InHalos[ihalo].Particles[ipart]=*it++;
  }
  assert(it==LocalParticles.end());
  vector <RemoteParticle_t>().swap(LocalParticles);
}

template <class Halo_T>
void ParticleExchanger_t<Halo_T>::Exchange()
{
/*  HBTInt ntot=0, ntot2=0;
#pragma omp parallel for reduction(+:ntot)
  for(HBTInt i=0;i<InHalos.size();i++)
    ntot+=InHalos[i].Particles.size();
*/
  while(true)
  {
    int npmax, np=CollectParticles();
//     ntot2+=np;
    MPI_Allreduce(&np, &npmax, 1, MPI_INT, MPI_MAX, world.Communicator);
    if(0==npmax) break; //loop till no particles to collect in all threads.
    SendParticles();
    QueryParticles();
    RecvParticles();
    RestoreParticles();
  }
//   assert(ntot==ntot2);

  for(auto &&h: InHalos)
    h.KickNullParticles();
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
void ParticleSnapshot_t::ExchangeHalos(MpiWorker_t& world, vector <Halo_T>& InHalos, vector<Halo_T>& OutHalos, MPI_Datatype MPI_Halo_Shell_Type, bool NeedParticleQuery) const
{
  typedef typename vector <Halo_T>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;

//   cout<<"Query particle..."<<flush;
  if(NeedParticleQuery)
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

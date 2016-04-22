#ifndef SNAPSHOT_H_INCLUDED
#define SNAPSHOT_H_INCLUDED

#include <iostream>
#include <sstream>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include <list>
#include <forward_list>

#include "datatypes.h"
#include "mymath.h"
#include "config_parser.h"
#include "snapshot_number.h"
#include "hash.h"
#include "mpi_wrapper.h"

#define NUMBER_OF_PARTICLE_TYPES 6
#define SNAPSHOT_HEADER_SIZE 256
struct IdRank_t
{
  HBTInt Id;
  int Rank;
  IdRank_t()=default;
  IdRank_t(HBTInt id, int rank): Id(id), Rank(rank)
  {
  }
};
#ifdef HBT_INT8
#define MPI_HBTRankPair MPI_LONG_INT
#else
#define MPI_HBTRankPair MPI_2INT
#endif
class Particle_t
{
public:
  HBTInt Id;
//   HBTInt ParticleIndex;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  HBTReal Mass;
  Particle_t()=default;
  Particle_t(HBTInt id): Id(id)
  {
  }
  void create_MPI_type(MPI_Datatype &MPI_HBTParticle_t);
};
extern ostream& operator << (ostream& o, Particle_t &p);
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

class SnapshotHeader_t
{
public:
  int      npart[NUMBER_OF_PARTICLE_TYPES];
  double   mass[NUMBER_OF_PARTICLE_TYPES];
  double   ScaleFactor;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned npartTotal[NUMBER_OF_PARTICLE_TYPES];  //differ from standard. to be able to hold large integers
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   OmegaM0;
  double   OmegaLambda0;
  double   HubbleParam; 
  char     fill[SNAPSHOT_HEADER_SIZE- NUMBER_OF_PARTICLE_TYPES*4- NUMBER_OF_PARTICLE_TYPES*8- 2*8- 2*4- NUMBER_OF_PARTICLE_TYPES*4- 2*4 - 4*8];  /* fills to 256 Bytes */
  void create_MPI_type(MPI_Datatype &dtype);
};

class Snapshot_t: public SnapshotNumber_t
{
public:
  /*epoch header*/
  HBTReal OmegaM0;
  HBTReal OmegaLambda0;
  HBTReal Hz; //current Hubble param in internal units
  HBTReal ScaleFactor;
  /*end epoch header*/
  
  Snapshot_t(): Hz(0.), ScaleFactor(0.), OmegaM0(0.), OmegaLambda0(0.), SnapshotNumber_t()
  {
  }
  Snapshot_t(const Snapshot_t & sn)
  {
	SetEpoch(sn);
  }
  void SetEpoch(HBTReal scalefactor, HBTReal omegaM0, HBTReal omegaLambda0)
  {
	ScaleFactor=scalefactor;
	OmegaM0=omegaM0;
	OmegaLambda0=omegaLambda0;
	Hz=PhysicalConst::H0 * sqrt(OmegaM0 / (ScaleFactor * ScaleFactor * ScaleFactor) 
	+ (1 - OmegaM0 - OmegaLambda0) / (ScaleFactor * ScaleFactor)
	+ OmegaLambda0);//Hubble param for the current catalogue;
  }
  void SetEpoch(const Snapshot_t & snap)
  {
	ScaleFactor=snap.ScaleFactor;
	Hz=snap.Hz;
	OmegaM0=snap.OmegaM0;
	OmegaLambda0=snap.OmegaLambda0;
  }
  virtual HBTInt size() const=0;
  virtual HBTInt GetId(HBTInt index) const
  {
	return index;
  }
  virtual const HBTxyz & GetComovingPosition(HBTInt index) const=0;
  virtual const HBTxyz & GetPhysicalVelocity(HBTInt index) const=0;
  virtual HBTReal GetMass(HBTInt index) const=0;
  void SphericalOverdensitySize(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <HBTReal> &RSorted, HBTReal ParticleMass) const;
  void SphericalOverdensitySize2(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <HBTReal> &RSorted, HBTReal ParticleMass) const;
  void HaloVirialFactors(HBTReal &virialF_tophat, HBTReal &virialF_b200, HBTReal &virialF_c200) const;
};
class SnapshotView_t: public Snapshot_t
{
public:
  HBTInt * Ids;
  HBTInt N;
  Snapshot_t & Snapshot;
  SnapshotView_t(vector <HBTInt> & ids, Snapshot_t & fullsnapshot): Ids(ids.data()), N(ids.size()), Snapshot(fullsnapshot)
  {
	SetEpoch(fullsnapshot);
  };
  SnapshotView_t(VectorView_t <HBTInt> &ids, Snapshot_t & fullsnapshot): Ids(ids.data()), N(ids.size()), Snapshot(fullsnapshot)
  {
	SetEpoch(fullsnapshot);
  };
  SnapshotView_t(HBTInt *ids, HBTInt n, Snapshot_t & fullsnapshot): Ids(ids), N(n), Snapshot(fullsnapshot)
  {
	SetEpoch(fullsnapshot);
  };
  void ReSize(HBTInt n)
  {
	N=n;
  }
  HBTInt size() const
  {
	return N;
  }
  HBTInt GetId(HBTInt i) const
  {
	return Snapshot.GetId(Ids[i]);
  }
  HBTReal GetMass(HBTInt i) const
  {
	return Snapshot.GetMass(Ids[i]);
  }
  const HBTxyz & GetPhysicalVelocity(HBTInt i) const
  {
	return Snapshot.GetPhysicalVelocity(Ids[i]);
  }
  const HBTxyz & GetComovingPosition(HBTInt i) const
  {
	return Snapshot.GetComovingPosition(Ids[i]);
  }
};

class ParticleSnapshot_t: public Snapshot_t
{
  typedef vector <HBTInt> IndexList_t;
  /*header extension*/
  bool NeedByteSwap;
  int IntTypeSize;
  int RealTypeSize;
  vector <HBTInt> NumberOfDMParticleInFiles;
  vector <HBTInt> OffsetOfDMParticleInFiles;
    
  HBTInt NumberOfParticles;
  FlatIndexTable_t<HBTInt, HBTInt> FlatHash;
  MappedIndexTable_t<HBTInt, HBTInt> MappedHash;
  IndexTable_t<HBTInt, HBTInt> *ParticleHash;
  
  void ReadFile(int ifile);
  void LoadHeader(int ifile=0);
  bool ReadFileHeader(FILE *fp, SnapshotHeader_t &header);
  HBTInt ReadNumberOfDMParticles(int ifile);
  size_t SkipBlock(FILE *fp);
  void ExchangeParticles(MpiWorker_t &world);
public:
  SnapshotHeader_t Header;
  vector <Particle_t> Particles;
  
  ParticleSnapshot_t(): Snapshot_t(), Header(), Particles(), NumberOfParticles(0), ParticleHash(), MappedHash(), FlatHash()
  {
	NeedByteSwap=false;
	IntTypeSize=0;
	RealTypeSize=0;
	if(HBTConfig.ParticleIdNeedHash)
	  ParticleHash=&MappedHash;
	else
	  ParticleHash=&FlatHash;
  }
  ~ParticleSnapshot_t()
  {
	Clear();//not necessary
  }
  void FillParticleHash();
  void ClearParticleHash();
  void GetFileName(int ifile, string &filename);
  
  HBTInt size() const;
  HBTInt GetId(HBTInt index) const;
  HBTInt GetIndex(HBTInt particle_id) const;
  HBTInt GetIndex(Particle_t & particle) const;
  const HBTxyz & GetComovingPosition(HBTInt index) const;
  const HBTxyz & GetPhysicalVelocity(HBTInt index) const;
  HBTReal GetMass(HBTInt index) const;
  
  void Load(MpiWorker_t &world, int snapshot_index, bool fill_particle_hash=true);
  void Clear();
  
  void AveragePosition(HBTxyz & CoM, const HBTInt Particles[], HBTInt NumPart) const; 
  void AverageVelocity(HBTxyz & CoV, const HBTInt Particles[], HBTInt NumPart) const;
  
  void MpiGetParticles(MpiWorker_t &world, vector <RemoteParticle_t> &particles) const;
  template <class Halo_T>
  void ExchangeHalos(MpiWorker_t &world, vector <Halo_T> & InHalos, vector <Halo_T> & OutHalos, MPI_Datatype MPI_Halo_Shell_Type) const;
};
inline HBTInt ParticleSnapshot_t::size() const
{
  return Particles.size();
}
inline HBTInt ParticleSnapshot_t::GetId(HBTInt index) const
{
  return Particles[index].Id;
}
inline HBTInt ParticleSnapshot_t::GetIndex(HBTInt particle_id) const
{
  return ParticleHash->GetIndex(particle_id);
}
inline HBTInt ParticleSnapshot_t::GetIndex(Particle_t & particle) const
{
  return ParticleHash->GetIndex(particle.Id);
}
inline const HBTxyz& ParticleSnapshot_t::GetComovingPosition(HBTInt index) const
{
  return Particles[index].ComovingPosition;
}
inline const HBTxyz& ParticleSnapshot_t::GetPhysicalVelocity(HBTInt index) const
{
  return Particles[index].PhysicalVelocity;
}
inline HBTReal ParticleSnapshot_t::GetMass(HBTInt index) const
{/*
  if(Header.mass[1])
	return Header.mass[1];
  else*/
	return Particles[index].Mass;
}

extern void AveragePosition(HBTxyz& CoM, const Particle_t Particles[], HBTInt NumPart);
extern void AverageVelocity(HBTxyz& CoV, const Particle_t Particles[], HBTInt NumPart);

extern void create_Mpi_RemoteParticleType(MPI_Datatype& dtype, bool IdOnly=false);
class ParticleExchanger_t
{
  int PrevRank, NextRank;
  int iloop, nloop;
  int iloop_debug, iloop_sent;
  const HBTInt EndParticleId;
  const int TagQuery;
  const int maxbuffersize;
  vector <RemoteParticle_t> sendbuffer, recvbuffer;
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
  ParticleExchanger_t(MpiWorker_t &_world, const ParticleSnapshot_t &_snap, HaloParticleIterator_t &particle_it): world(_world), snap(_snap), iloop(0), EndParticleId(-1), TagQuery(1), maxbuffersize(100), ReqSend(MPI_REQUEST_NULL), ChannelIsClean(1), iloop_debug(1), iloop_sent(0), sendbuffer(maxbuffersize), recvbuffer(maxbuffersize)
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

inline bool CompareRank(const IdRank_t &a, const IdRank_t &b)
{
  return (a.Rank<b.Rank);
}

template <class Halo_T>
void ParticleSnapshot_t::ExchangeHalos(MpiWorker_t& world, vector <Halo_T>& InHalos, vector<Halo_T>& OutHalos, MPI_Datatype MPI_Halo_Shell_Type) const
{
  typedef typename vector <Halo_T>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;
  
  cout<<"Query particle..."<<flush;
  vector <IdRank_t>TargetRank;
  {//query particles
	ParticleIterator_t InParticles(InHalos.begin(), InHalos.end());
	ParticleExchanger_t Exchanger(world, *this, InParticles);
	Exchanger.Exchange();
	Exchanger.UnPackHaloParticles(InHalos, TargetRank);
  }
  
  cout<<"sending shells...\n";
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
	
	cout<<"sending particles...";
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
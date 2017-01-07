#include <iostream>
#include <new>
#include <algorithm>
#include <unordered_set>
#include <omp.h>

#include "datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"

void Subhalo_t::UpdateTrack(const Snapshot_t &epoch)
{
  if(TrackId==SpecialConst::NullTrackId) return;
  
  if(0==Rank) SnapshotIndexOfLastIsolation=epoch.GetSnapshotIndex();
  if(Mbound>=LastMaxMass) 
  {
	SnapshotIndexOfLastMaxMass=epoch.GetSnapshotIndex();
	LastMaxMass=Mbound;
  }
}
HBTReal Subhalo_t::KineticDistance(const Halo_t &halo, const Snapshot_t &epoch)
{
  HBTxyz dv;
  epoch.RelativeVelocity(ComovingAveragePosition, PhysicalAverageVelocity, halo.ComovingAveragePosition, halo.PhysicalAverageVelocity, dv);
  return VecNorm(dv);
}
void MemberShipTable_t::ResizeAllMembers(size_t n)
{
  auto olddata=AllMembers.data();
  AllMembers.resize(n);
  size_t offset=AllMembers.data()-olddata;
  if(offset)
  {
	for(HBTInt i=0;i<Mem_SubGroups.size();i++)
	  Mem_SubGroups[i].Bind(Mem_SubGroups[i].data()+offset);
  }
}
void MemberShipTable_t::Init(const HBTInt nhalos, const HBTInt nsubhalos, const float alloc_factor)
{
  Mem_SubGroups.clear();
  Mem_SubGroups.resize(nhalos+1);
  SubGroups.Bind(nhalos, Mem_SubGroups.data()+1);

  AllMembers.clear();
  AllMembers.reserve(nsubhalos*alloc_factor); //allocate more for seed haloes.
  AllMembers.resize(nsubhalos);
}
void MemberShipTable_t::BindMemberLists()
{
  HBTInt offset=0;
  for(HBTInt i=0;i<Mem_SubGroups.size();i++)
  {
	Mem_SubGroups[i].Bind(Mem_SubGroups[i].size(), &(AllMembers[offset]));
	offset+=Mem_SubGroups[i].size();
	Mem_SubGroups[i].ReBind(0);
  }
}
void MemberShipTable_t::CountMembers(const SubhaloList_t& Subhalos, bool include_orphans)
{//todo: parallelize this..
  if(include_orphans)
  {
	for(HBTInt subid=0;subid<Subhalos.size();subid++)
	  SubGroups[Subhalos[subid].HostHaloId].IncrementBind();
  }
  else
  {
	for(HBTInt subid=0;subid<Subhalos.size();subid++)
	  if(Subhalos[subid].Nbound>1)
		SubGroups[Subhalos[subid].HostHaloId].IncrementBind();
  }
}
void MemberShipTable_t::FillMemberLists(const SubhaloList_t& Subhalos, bool include_orphans)
{//fill with local subhaloid
  if(include_orphans)
  {
	for(HBTInt subid=0;subid<Subhalos.size();subid++)
	  SubGroups[Subhalos[subid].HostHaloId].PushBack(subid);
  }
  else
  {
	for(HBTInt subid=0;subid<Subhalos.size();subid++)
	  if(Subhalos[subid].Nbound>1)
		SubGroups[Subhalos[subid].HostHaloId].PushBack(subid);
  }
}
struct CompareMass_t
{
  const SubhaloList_t * Subhalos;
  CompareMass_t(const SubhaloList_t & subhalos)
  {
	Subhalos=&subhalos;
  }
  bool operator () (const HBTInt & i, const HBTInt & j)
  {
	return (*Subhalos)[i].Mbound>(*Subhalos)[j].Mbound;
  }
};
void MemberShipTable_t::SortMemberLists(const SubhaloList_t & Subhalos)
{
  CompareMass_t compare_mass(Subhalos);
#pragma omp for
  for(HBTInt i=-1;i<SubGroups.size();i++)
	std::sort(SubGroups[i].begin(), SubGroups[i].end(), compare_mass);
}
void MemberShipTable_t::SortSatellites(const SubhaloList_t & Subhalos)
/*central subhalo not changed*/
{
  CompareMass_t compare_mass(Subhalos);
  for(HBTInt i=0;i<SubGroups.size();i++)
	std::sort(SubGroups[i].begin()+1, SubGroups[i].end(), compare_mass);
}
void MemberShipTable_t::AssignRanks(SubhaloList_t& Subhalos)
{
  //field subhaloes
  {
	MemberList_t & SubGroup=SubGroups[-1];
  #pragma omp for
	for(HBTInt i=0;i<SubGroup.size();i++)
	  Subhalos[SubGroup[i]].Rank=0;
  }
#pragma omp for
  for(HBTInt haloid=0;haloid<SubGroups.size();haloid++)
  {
	MemberList_t & SubGroup=SubGroups[haloid];
	for(HBTInt i=0;i<SubGroup.size();i++)
	  Subhalos[SubGroup[i]].Rank=i;
  }
}
void MemberShipTable_t::CountEmptyGroups()
{
static HBTInt nbirth;
#pragma omp single
nbirth=0;
#pragma omp for reduction(+:nbirth)
  for(HBTInt hostid=0;hostid<SubGroups.size();hostid++)
	if(SubGroups[hostid].size()==0)  nbirth++;
#pragma omp single
NBirth=nbirth;	
}
/*
inline bool SubhaloSnapshot_t::CompareHostAndMass(const HBTInt& subid_a, const HBTInt& subid_b)
{//ascending in host id, descending in mass inside each host, and put NullHaloId to the beginning.
  Subhalo_t a=Subhalos[subid_a], b=Subhalos[subid_b];
  
  if(a.HostHaloId==b.HostHaloId) return (a.Nbound>b.Nbound);
  
  return (a.HostHaloId<b.HostHaloId); //(a.HostHaloId!=SpecialConst::NullHaloId)&&
}*/
void MemberShipTable_t::Build(const HBTInt nhalos, const SubhaloList_t & Subhalos, bool include_orphans)
{
  #pragma omp single
  {
  Init(nhalos, Subhalos.size());
  CountMembers(Subhalos, include_orphans);
  BindMemberLists();
  FillMemberLists(Subhalos, include_orphans);
  }
  SortMemberLists(Subhalos);
  CountEmptyGroups();
//   std::sort(AllMembers.begin(), AllMembers.end(), CompareHostAndMass);
}

inline HBTInt GetLocalHostId(HBTInt pid, const HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap)
{
	HBTInt hostid=halo_snap.ParticleHash.GetIndex(pid);
	if(hostid<0)//not in the haloes, =-1
	{
	  if(part_snap.GetIndex(pid)==SpecialConst::NullParticleId)
		hostid--;//not in this snapshot either, =-2
	}
	return hostid;
}
void FindLocalHosts(const HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap, vector <Subhalo_t> & Subhalos, vector <Subhalo_t> &LocalSubhalos)
{
  #pragma omp parallel for 
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	if(Subhalos[subid].Particles.size())
	  Subhalos[subid].HostHaloId=GetLocalHostId(Subhalos[subid].Particles[0].Id, halo_snap, part_snap);
	else
	  Subhalos[subid].HostHaloId=-1;
  }
  
  HBTInt nsub=0;
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	if(Subhalos[subid].HostHaloId<0&&Subhalos[subid].Particles.size())//only move nonempty subhalos
	{
	  if(subid>nsub)
		Subhalos[nsub]=move(Subhalos[subid]);//there should be a default move assignement operator.
	  nsub++;
	}
	else
	  LocalSubhalos.push_back(move(Subhalos[subid]));
  }
  Subhalos.resize(nsub);
}

void FindOtherHosts(MpiWorker_t &world, int root, const HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap, VectorView_t <Subhalo_t> &Subhalos, vector <Subhalo_t> &LocalSubhalos, MPI_Datatype MPI_Subhalo_Shell_Type)
/*scatter Subhalos from process root to LocalSubhalos in every other process
 Note Subalos are "moved", so are in a unspecified state upon return.*/
{
  int thisrank=world.rank();
  vector <HBTInt> TrackParticleIds;
  HBTInt NumSubhalos;
  
  //broadcast trackparticles
  if(thisrank==root) 
  {
	NumSubhalos=Subhalos.size();
	if(NumSubhalos>INT_MAX)
	  throw runtime_error("Error: in FindOtherHosts(), sending more subhaloes than INT_MAX will cause MPI message to overflow. Please try more MPI threads. aborting.\n");
  }
  MPI_Bcast(&NumSubhalos, 1, MPI_HBT_INT, root, world.Communicator);
  TrackParticleIds.resize(NumSubhalos);
  if(thisrank==root)
  {
	for(HBTInt i=0;i<Subhalos.size();i++)
	  TrackParticleIds[i]=Subhalos[i].Particles[0].Id;
  }
  MPI_Bcast(TrackParticleIds.data(), NumSubhalos, MPI_HBT_INT, root, world.Communicator);
  
  //find hosts
  vector <IdRank_t> LocalHostIds(NumSubhalos), GlobalHostIds(NumSubhalos);
  if(thisrank==root)
  {
	#pragma omp parallel for if(NumSubhalos>20)
	for(HBTInt i=0;i<NumSubhalos;i++)
	{
	  LocalHostIds[i].Id=Subhalos[i].HostHaloId;//already found previously
	  LocalHostIds[i].Rank=thisrank;
	}
  }
  else
  {
	#pragma omp parallel for if(NumSubhalos>20)
	for(HBTInt i=0;i<NumSubhalos;i++)
	{
	  LocalHostIds[i].Id=GetLocalHostId(TrackParticleIds[i], halo_snap, part_snap);
	  LocalHostIds[i].Rank=thisrank;
	}
  }

  MPI_Allreduce(LocalHostIds.data(), GlobalHostIds.data(), NumSubhalos, MPI_HBTRankPair, MPI_MAXLOC, world.Communicator);

  //scatter free subhaloes from root to everywhere
  //send particles and nests; no scatterw, do it manually
  MPI_Datatype MPI_HBT_Particle;
  Particle_t().create_MPI_type(MPI_HBT_Particle);
  vector <vector<int> > SendSizes(world.size()), SendNestSizes(world.size());
  vector <MPI_Request> Req0,Req1, ReqNest0, ReqNest1;
  if(thisrank==root)
  {
	vector <vector <MPI_Aint> > SendBuffers(world.size()), SendNestBuffers(world.size());
	for(HBTInt subid=0;subid<NumSubhalos;subid++)//packing
	{
	  int rank=GlobalHostIds[subid].Rank;
	  auto & Particles=Subhalos[subid].Particles;
	  MPI_Aint p;
	  MPI_Address(Particles.data(),&p);
	  SendBuffers[rank].push_back(p);
	  SendSizes[rank].push_back(Particles.size());
	  
	  auto & Nest=Subhalos[subid].NestedSubhalos;
	  MPI_Address(Nest.data(),&p);
	  SendNestBuffers[rank].push_back(p);
	  SendNestSizes[rank].push_back(Nest.size());
	}
	Req0.resize(world.size());
	Req1.resize(world.size());
	ReqNest0.resize(world.size());
	ReqNest1.resize(world.size());
	for(int rank=0;rank<world.size();rank++)
	{
	  {
	  MPI_Isend(SendSizes[rank].data(), SendSizes[rank].size(), MPI_INT, rank, 0, world.Communicator, &Req0[rank]);
	  MPI_Datatype SendType;
	  MPI_Type_create_hindexed(SendSizes[rank].size(), SendSizes[rank].data(), SendBuffers[rank].data(), MPI_HBT_Particle, &SendType);
	  MPI_Type_commit(&SendType);
	  MPI_Isend(MPI_BOTTOM, 1, SendType, rank, 1, world.Communicator, &Req1[rank]);
	  MPI_Type_free(&SendType);
	  }
	  {
	  MPI_Isend(SendNestSizes[rank].data(), SendNestSizes[rank].size(), MPI_INT, rank, 0, world.Communicator, &ReqNest0[rank]);
	  MPI_Datatype SendNestType;
	  MPI_Type_create_hindexed(SendNestSizes[rank].size(), SendNestSizes[rank].data(), SendNestBuffers[rank].data(), MPI_HBT_INT, &SendNestType);
	  MPI_Type_commit(&SendNestType);
	  MPI_Isend(MPI_BOTTOM, 1, SendNestType, rank, 1, world.Communicator, &ReqNest1[rank]);
	  MPI_Type_free(&SendNestType);
	  }
	}
  }
  //receive on every process, including root
	vector <MPI_Aint> ReceiveBuffer, ReceiveNestBuffer;
	vector <int> ReceiveSize, ReceiveNestSize;
	int NumNewSubs;
	MPI_Status stat;
	MPI_Probe(root, 0, world.Communicator, &stat);
	MPI_Get_count(&stat, MPI_INT, &NumNewSubs);
	ReceiveSize.resize(NumNewSubs);
	MPI_Recv(ReceiveSize.data(), NumNewSubs, MPI_INT, root, 0, world.Communicator, &stat);
	LocalSubhalos.resize(LocalSubhalos.size()+NumNewSubs);
	auto NewSubhalos=LocalSubhalos.end()-NumNewSubs;
	ReceiveBuffer.resize(NumNewSubs);
	for(int i=0;i<NumNewSubs;i++)
	{
	  auto &Particles=NewSubhalos[i].Particles;
	  Particles.resize(ReceiveSize[i]);
	  MPI_Aint p;
	  MPI_Address(Particles.data(),&p);
	  ReceiveBuffer[i]=p;
	}
	MPI_Datatype ReceiveType;
	MPI_Type_create_hindexed(NumNewSubs, ReceiveSize.data(), ReceiveBuffer.data(), MPI_HBT_Particle, &ReceiveType);
	MPI_Type_commit(&ReceiveType);
	MPI_Recv(MPI_BOTTOM, 1, ReceiveType, root, 1, world.Communicator, &stat);
	MPI_Type_free(&ReceiveType);
	
	MPI_Type_free(&MPI_HBT_Particle);
	
	ReceiveNestSize.resize(NumNewSubs);
	MPI_Recv(ReceiveNestSize.data(), NumNewSubs, MPI_INT, root, 0, world.Communicator, &stat);
	ReceiveNestBuffer.resize(NumNewSubs);
	for(int i=0;i<NumNewSubs;i++)
	{
	  auto &Nest=NewSubhalos[i].NestedSubhalos;
	  Nest.resize(ReceiveNestSize[i]);
	  MPI_Aint p;
	  MPI_Address(Nest.data(),&p);
	  ReceiveNestBuffer[i]=p;
	}
	MPI_Datatype ReceiveNestType;
	MPI_Type_create_hindexed(NumNewSubs, ReceiveNestSize.data(), ReceiveNestBuffer.data(), MPI_HBT_INT, &ReceiveNestType);
	MPI_Type_commit(&ReceiveNestType);
	MPI_Recv(MPI_BOTTOM, 1, ReceiveNestType, root, 1, world.Communicator, &stat);
	MPI_Type_free(&ReceiveNestType);
	
	if(thisrank==root)
	{
// 	  vector <MPI_Status> stats(world.size());
	  MPI_Waitall(world.size(), Req0.data(), MPI_STATUS_IGNORE);
	  MPI_Waitall(world.size(), Req1.data(), MPI_STATUS_IGNORE);
	  MPI_Waitall(world.size(), ReqNest0.data(), MPI_STATUS_IGNORE);
	  MPI_Waitall(world.size(), ReqNest1.data(), MPI_STATUS_IGNORE);
	}
	
  //copy other properties
    vector <int> Counts(world.size()), Disps(world.size());
	vector <Subhalo_t> TmpHalos;
	if(world.rank()==root)
	{//reuse GlobalHostIds for sorting
	  for(HBTInt subid=0;subid<Subhalos.size();subid++)
	  {
		Subhalos[subid].HostHaloId=GlobalHostIds[subid].Id;
// 		assert(GlobalHostIds[subid].n>=-1);
		GlobalHostIds[subid].Id=subid;
	  }
	  stable_sort(GlobalHostIds.begin(), GlobalHostIds.end(), CompareRank);
	  TmpHalos.resize(Subhalos.size());
	  for(HBTInt subid=0;subid<Subhalos.size();subid++)
		TmpHalos[subid]=move(Subhalos[GlobalHostIds[subid].Id]);
// 		Subhalos[GlobalHostIds[subid].n].MoveTo(TmpHalos[subid], false);
	  for(int rank=0;rank<world.size();rank++)
		Counts[rank]=SendSizes[rank].size();
	  CompileOffsets(Counts, Disps);
	}
	MPI_Scatterv(TmpHalos.data(), Counts.data(), Disps.data(), MPI_Subhalo_Shell_Type, &NewSubhalos[0], NumNewSubs, MPI_Subhalo_Shell_Type, root, world.Communicator);
}
void FindOtherHostsSafely(MpiWorker_t &world, int root, const HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap, vector <Subhalo_t> &Subhalos, vector <Subhalo_t> &LocalSubhalos, MPI_Datatype MPI_Subhalo_Shell_Type)
/*break Subhalos into small chunks and then FindOtherHosts() for them, to avoid overflow in MPI message size*/
{
  const int MaxChunkSize=1024*1024;
  int flagstop=0;
  VectorView_t <Subhalo_t> HaloChunk;
  if(world.rank()==root)
  {
	HBTInt offset=0, chunksize=0;
	for(HBTInt i=0;i<Subhalos.size();i++)
	{
	  chunksize+=Subhalos[i].Particles.size();
	  if(chunksize>=MaxChunkSize)
	  {//send buffer
		HaloChunk.Bind(i-offset, &Subhalos[offset]);
		if(HaloChunk.size())
		{
		  MPI_Bcast(&flagstop, 1, MPI_INT, root, world.Communicator);
		  FindOtherHosts(world, root, halo_snap, part_snap, HaloChunk, LocalSubhalos, MPI_Subhalo_Shell_Type);
		}
		//reset buffer
		offset=i;
		chunksize=Subhalos[i].Particles.size();//the current halo
	  }
	}
	//remaining ones
	if(Subhalos.size()>offset)
	{
	  HaloChunk.Bind(Subhalos.size()-offset, &Subhalos[offset]);
	  MPI_Bcast(&flagstop, 1, MPI_INT, root, world.Communicator);
	  FindOtherHosts(world, root, halo_snap, part_snap, HaloChunk, LocalSubhalos, MPI_Subhalo_Shell_Type);
	}
	flagstop=1;
	MPI_Bcast(&flagstop, 1, MPI_INT, root, world.Communicator);
  }
  else
  {
	while(true)
	{
	  MPI_Bcast(&flagstop, 1, MPI_INT, root, world.Communicator);
	  if(flagstop) 
		break;
	  else
	    FindOtherHosts(world, root, halo_snap, part_snap, HaloChunk, LocalSubhalos, MPI_Subhalo_Shell_Type);
	}
  }
}
void SubhaloSnapshot_t::AssignHosts(MpiWorker_t &world, HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap)
/* find host haloes for subhaloes, and build MemberTable. Each subhalo is moved to the processor of its host halo, with its HostHaloId set to the local haloid of the host*/
{
  ParallelizeHaloes=halo_snap.NumPartOfLargestHalo<0.1*halo_snap.TotNumberOfParticles;//no dominating objects
  
  //exchange subhaloes according to hosts
  vector <Subhalo_t> LocalSubhalos;
  LocalSubhalos.reserve(Subhalos.size());
  halo_snap.FillParticleHash();
  FindLocalHosts(halo_snap, part_snap, Subhalos, LocalSubhalos);
  for(int rank=0;rank<world.size();rank++)
	FindOtherHostsSafely(world, rank, halo_snap, part_snap, Subhalos, LocalSubhalos, MPI_HBT_SubhaloShell_t);
  Subhalos.swap(LocalSubhalos);
  halo_snap.ClearParticleHash();
  
  MemberTable.Build(halo_snap.Halos.size(), Subhalos, true);
//   MemberTable.AssignRanks(Subhalos); //not needed here
}

void SubhaloSnapshot_t::DecideCentrals(const HaloSnapshot_t &halo_snap)
/* to select central subhalo according to KineticDistance, and move each central to the beginning of each list in MemberTable*/
{
#pragma omp for
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
  {
	MemberShipTable_t::MemberList_t &List=MemberTable.SubGroups[hostid];
	if(List.size()>1)
	{	  
	  int n_major;
	  HBTInt MassLimit=Subhalos[List[0]].Mbound*HBTConfig.MajorProgenitorMassRatio;
	  for(n_major=1;n_major<List.size();n_major++)
		if(Subhalos[List[n_major]].Mbound<MassLimit) break;
	  if(n_major>1)
	  {
		HBTReal dmin=Subhalos[List[0]].KineticDistance(halo_snap.Halos[hostid], *this);
		int icenter=0;
		for(int i=1;i<n_major;i++)
		{
		  HBTReal d=Subhalos[List[i]].KineticDistance(halo_snap.Halos[hostid], *this);
		  if(dmin>d)
		  {
			dmin=d;
			icenter=i;
		  }
		}
		if(icenter)  swap(List[0], List[icenter]);
	  }
	}
  }
}

void SubhaloSnapshot_t::FeedCentrals(HaloSnapshot_t& halo_snap)
/* replace centrals with host particles;
 * create a new central if there is none
 * initialize new central with host halo center coordinates
 * halo_snap is rendered unspecified upon return (its particles have been swapped to subhaloes).
 */
{
  static HBTInt Npro;
  #pragma omp single
  {
  Npro=Subhalos.size();
  //Subhalos.reserve(Snapshot->size()*0.1);//reserve enough	branches.......
  Subhalos.resize(Npro+MemberTable.NBirth);
  }
  #pragma omp for
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
  { 
	MemberShipTable_t::MemberList_t & Members=MemberTable.SubGroups[hostid];
	auto &Host=halo_snap.Halos[hostid];
	if(0==Members.size()) //create a new sub
	{
	  HBTInt subid;
	  #pragma omp critical(AddNewSub) //maybe consider ordered for here..
	  {
		subid=Npro++;
	  }
	  auto &central=Subhalos[subid];
	  central.HostHaloId=hostid;
	  copyHBTxyz(central.ComovingAveragePosition, Host.ComovingAveragePosition); 
	  copyHBTxyz(central.PhysicalAverageVelocity, Host.PhysicalAverageVelocity);
	  central.Particles.swap(Host.Particles);
	  central.Nbound=central.Particles.size();//init Nbound to source size.
	  central.SnapshotIndexOfBirth=SnapshotIndex;
	}
	else
	{
	  auto &central=Subhalos[Members[0]];
	  assert(central.Particles.size());
	  central.Particles.swap(Host.Particles); //reuse the halo particles
	  central.Nbound=central.Particles.size();
	  {
	    auto mostbndid=Host.Particles[0].Id; 
	    for(auto & p: central.Particles)
	    if(p.Id==mostbndid)//swap previous mostbound particle to the beginning
	    {
		  swap(p, central.Particles[0]);
		  break;
	    }
	  }
	}
  }
//   #pragma omp single
//   halo_snap.Clear();//to avoid misuse
}
void SubhaloSnapshot_t::PrepareCentrals(MpiWorker_t &world, HaloSnapshot_t &halo_snap)
{
  #pragma omp parallel
  {
  DecideCentrals(halo_snap);
  FeedCentrals(halo_snap);
  }
  NestSubhalos(world);
#ifndef INCLUSIVE_MASS
  MaskSubhalos();
#endif
}

void SubhaloSnapshot_t::RegisterNewTracks(MpiWorker_t &world)
/*assign trackId to new bound ones, remove unbound ones, and rebuild membership*/
{
  HBTInt NumSubMax=Subhalos.size(), NumSubOld=NumSubMax-MemberTable.NBirth;
  HBTInt NumSubNew=NumSubOld;
  for(HBTInt i=NumSubNew;i<NumSubMax;i++)
  {
	if(Subhalos[i].Nbound>1)
	{
	  if(i!=NumSubNew)
		Subhalos[NumSubNew]=move(Subhalos[i]);
	  NumSubNew++;
	}
  }
  Subhalos.resize(NumSubNew);
  MemberTable.Build(MemberTable.SubGroups.size(), Subhalos, true);//rebuild membership with new subs and also include orphans this time.
  MemberTable.NFake=NumSubMax-NumSubNew;
  MemberTable.NBirth=NumSubNew-NumSubOld;
  
  //now assign a global TrackId
  HBTInt TrackIdOffset, NBirth=MemberTable.NBirth, GlobalNumberOfSubs;
  MPI_Allreduce(&NumSubOld, &GlobalNumberOfSubs, 1, MPI_HBT_INT, MPI_SUM, world.Communicator); 
  MPI_Scan(&NBirth, &TrackIdOffset, 1, MPI_HBT_INT, MPI_SUM, world.Communicator); 
  TrackIdOffset=TrackIdOffset+GlobalNumberOfSubs-NBirth;
  for(HBTInt i=NumSubOld;i<NumSubNew;i++)
	Subhalos[i].TrackId=TrackIdOffset++;
}
void SubhaloSnapshot_t::PurgeMostBoundParticles()
/* fix the possible issue that the most-bound particle of a subhalo might belong to a sub-sub.
 * this is achieved by masking particles from smaller subs and promote the remaining most-bound particle in each subhalo.
 * orphan galaxies are not considered.
 * the shift in center position does not affect the total angular momentum.
 */
{
  #pragma omp for
  for(HBTInt i=-1;i<MemberTable.SubGroups.size();i++)
  {
	auto &Group=MemberTable.SubGroups[i];
	unordered_set <HBTInt> ExclusionList;
	{
	HBTInt np=0;
	for(auto &&subid: Group)
	  if(Subhalos[subid].Nbound>1) np+=Subhalos[subid].Nbound;
	ExclusionList.reserve(np);
	}
	for(HBTInt j=Group.size()-1;j>=0;j--)
	{
	  auto & subhalo=Subhalos[Group[j]];
	  if(subhalo.Nbound>1)
	  {
		for(auto & p: subhalo.Particles)
		{
		  if(ExclusionList.find(p.Id)==ExclusionList.end())
		  {
			if(&p!=&subhalo.Particles[0])
			{
			  copyHBTxyz(subhalo.ComovingMostBoundPosition, p.ComovingPosition);
			  copyHBTxyz(subhalo.PhysicalMostBoundVelocity, p.PhysicalVelocity);
			  swap(subhalo.Particles[0], p);
			  break;
			}
		  }
		}
		for(auto && p: subhalo.Particles)
		  ExclusionList.insert(p.Id);//alternative: only filter most-bounds
	  }
	}
  }
}

void SubhaloSnapshot_t::LocalizeNestedIds(MpiWorker_t &world)
/*convert TrackIds of NestedSubhalos to local index, and move non-local nestedsubhalos to their host processors as DissociatedTrack*/
{
  TrackKeyList_t Ids(*this); 
  MappedIndexTable_t<HBTInt, HBTInt> TrackHash;
  TrackHash.Fill(Ids, SpecialConst::NullTrackId);
  
  //collect lost tracks
  vector <HBTInt> DissociatedTracks;
  for(auto &&subhalo: Subhalos)
  {
    auto & nests=subhalo.NestedSubhalos;
    auto it_begin=nests.begin();
    auto it_save=it_begin, it=it_begin;
    for(;it!=nests.end();++it)
    {
      HBTInt subid=TrackHash.GetIndex(*it);
      if(subid==SpecialConst::NullTrackId)
      {
	DissociatedTracks.push_back(*it);
      }
      else
      {
	*it_save=subid;
	++it_save;
      }
    }
    nests.resize(it_save-it_begin);
  }
  
  //distribute, locate and levelup DissociatedTracks
  vector <HBTInt> ReceivedTracks;
  for(int root=0;root<world.size();root++)
  {
    HBTInt stacksize=DissociatedTracks.size();
    MPI_Bcast(&stacksize, 1, MPI_HBT_INT, root, world.Communicator);
    ReceivedTracks.resize(stacksize);
    MyBcast<HBTInt, vector <HBTInt>::iterator, vector <HBTInt>::iterator>(world, DissociatedTracks.begin(), ReceivedTracks.begin(), stacksize, MPI_HBT_INT, root);
    if(world.rank()!=root)
    {
      for(auto & tid: ReceivedTracks)
      {
	auto subid=TrackHash.GetIndex(tid);
	if(subid!=SpecialConst::NullTrackId)//located
	  Subhalos[subid].Rank=0; //level up this DissociatedTrack
      }
    }
  }
}
void SubhaloSnapshot_t::GlobalizeNestedIds()
/*translate subhalo index to trackIds for nestedsubhalos*/
{
    for(auto &subhalo: Subhalos)
      for(auto &&subid: subhalo.NestedSubhalos)
	subid=Subhalos[subid].TrackId;
}
void SubhaloSnapshot_t::NestSubhalos(MpiWorker_t &world)
{
  LocalizeNestedIds(world);
  LevelUpDetachedSubhalos();
  //collect detached(head) subhalos
  #pragma omp single
  MemberTable.SubGroupsOfHeads.clear();
  MemberTable.SubGroupsOfHeads.resize(MemberTable.SubGroups.size());
  #pragma omp parallel for
  for(HBTInt haloid=0;haloid<MemberTable.SubGroups.size();haloid++)
  {
    auto &subgroup=MemberTable.SubGroups[haloid];
    for(HBTInt i=0;i<subgroup.size();i++)
    {
      auto subid=subgroup[i];
      if(Subhalos[subid].Rank==0)
	MemberTable.SubGroupsOfHeads[haloid].push_back(subid);
    }
  }  
}

void SubhaloSnapshot_t::ExtendCentralNest()
{
  #pragma omp for
  for(HBTInt haloid=0;haloid<MemberTable.SubGroups.size();haloid++)
  {
    auto &subgroup=MemberTable.SubGroups[haloid];
    auto &heads=MemberTable.SubGroupsOfHeads[haloid];
    if(subgroup.size()<=1) continue;
    auto &central=Subhalos[subgroup[0]];
    if(central.Rank==0)//already a central
    {
      central.NestedSubhalos.reserve(central.NestedSubhalos.size()+heads.size()-1);
      for(auto &&h: heads)
	if(h!=subgroup[0]) central.NestedSubhalos.push_back(h);
    }
    else
    {
      central.Rank=0; //promote to central
      for(auto &&h: heads)
	Subhalos[h].LevelUpDetachedMembers(Subhalos);
      central.NestedSubhalos.insert(central.NestedSubhalos.end(), heads.begin(), heads.end());
    }
  }
#pragma omp single
  MemberTable.SubGroupsOfHeads.clear();
}

void SubhaloSnapshot_t::LevelUpDetachedSubhalos()
/*
 * assign rank=0 to subhaloes that has drifted away from the hosthalo of its host-subhalo.
 */
{
  vector <char> IsHeadSub(Subhalos.size());
//record head list first, since the ranks are modified during LevelUpDetachedMembers().
#pragma omp parallel
  {
  #pragma omp for
    for(HBTInt subid=0; subid<Subhalos.size(); subid++)
	IsHeadSub[subid]=(Subhalos[subid].Rank==0);
    
  //promote centrals to detached
  #pragma omp for
    for(HBTInt haloid=0;haloid<MemberTable.SubGroups.size();haloid++)
    {
      auto &subgroup=MemberTable.SubGroups[haloid];
      if(subgroup.size())
	Subhalos[subgroup[0]].Rank=0;
    }
    {
    auto &subgroup=MemberTable.SubGroups[-1];
  #pragma omp for
    for(HBTInt i=0; i<subgroup.size();i++)//break up all field subhalos
      Subhalos[subgroup[i]].Rank=0;
    }
    //TODO: break up all orphans as well?
    
  #pragma omp for
    for(HBTInt subid=0;subid<Subhalos.size();subid++)
      if(IsHeadSub[subid])
	Subhalos[subid].LevelUpDetachedMembers(Subhalos);
  }
}

void Subhalo_t::LevelUpDetachedMembers(vector <Subhalo_t> &Subhalos)
{
  HBTInt isave=0;
  for(HBTInt i=0;i<NestedSubhalos.size();i++)
  {
    auto subid=NestedSubhalos[i];
    if(Subhalos[subid].HostHaloId!=HostHaloId||Subhalos[subid].Rank==0)
    {
      if(Subhalos[subid].Rank)
	Subhalos[subid].Rank=0;
    }
    else
    {
      if(isave!=i)
	NestedSubhalos[isave]=subid;
      isave++;
    }
    Subhalos[subid].LevelUpDetachedMembers(Subhalos);//recursively level up members. Note this can be further improved: if its members didn't follow it but stayed in the original host, then they should be added to the current NestedSubhalos, instead of being leveled up to rank 0. probably not necessary, since you did not actually check the host-sub but only adopted the historical relation.
  }
  NestedSubhalos.resize(isave);//remove detached ones from the list
}

class SubhaloMasker_t
{
  unordered_set <HBTInt> ExclusionList;
public:
  SubhaloMasker_t(HBTInt np_guess)
  {
    ExclusionList.reserve(np_guess);
  }
  void Mask(HBTInt subid, vector <Subhalo_t> &Subhalos)
  {
    auto &subhalo=Subhalos[subid];
    for(auto nestedid: subhalo.NestedSubhalos)//TODO: do we have to do it recursively? satellites are already masked among themselves?
      Mask(nestedid, Subhalos);
	
	if(subhalo.Nbound<=1) return; //skip orphans
    
    auto it_begin=subhalo.Particles.begin(), it_save=it_begin;
    for(auto it=it_begin;it!=subhalo.Particles.end();++it)
    {
      auto insert_status=ExclusionList.insert(it->Id);
      if(insert_status.second)//inserted, meaning not excluded
      {
	if(it!=it_save)
	  *it_save=move(*it);
	++it_save;
      }
    }
    subhalo.Particles.resize(it_save-it_begin);
  }
};

void SubhaloSnapshot_t::MaskSubhalos()
{
  #pragma omp parallel for
  for(HBTInt i=0;i<MemberTable.SubGroups.size();i++)
  {
	auto &Group=MemberTable.SubGroups[i];
	if(Group.size()==0) continue;
	auto &central=Subhalos[Group[0]];
	auto &nest=central.NestedSubhalos;
	auto old_membercount=nest.size();
	auto &heads=MemberTable.SubGroupsOfHeads[i];
	//update central member list (append other heads except itself)
	nest.insert(nest.end(), heads.begin()+1, heads.end());
	SubhaloMasker_t Masker(central.Particles.size()*1.2);
	Masker.Mask(Group[0], Subhalos);
	nest.resize(old_membercount);//TODO: better way to do this? or do not change the nest for central?
  }
}

void SubhaloSnapshot_t::UpdateTracks(MpiWorker_t &world, const HaloSnapshot_t &halo_snap)
{
  /*renew ranks after unbinding*/
  RegisterNewTracks(world);//performance bottleneck here. no. just poor synchronization.
  #pragma omp parallel
  {
  MemberTable.SortMemberLists(Subhalos);//reorder, so the central might change if necessary
  ExtendCentralNest();
  MemberTable.AssignRanks(Subhalos);
#ifdef INCLUSIVE_MASS
  PurgeMostBoundParticles();
#endif
#pragma omp for
  for(HBTInt i=0;i<Subhalos.size();i++)
  {
	Subhalos[i].UpdateTrack(*this);
	HBTInt HostId=Subhalos[i].HostHaloId;
	if(HostId<0)
	  Subhalos[i].HostHaloId=-1;
	else
	  Subhalos[i].HostHaloId=halo_snap.Halos[HostId].HaloId;//restore global haloid
  }
  }
  GlobalizeNestedIds();
  #pragma omp parallel for if(ParallelizeHaloes)
  for(HBTInt i=0;i<Subhalos.size();i++)
  {
	Subhalos[i].CalculateProfileProperties(*this);
	Subhalos[i].CalculateShape();
  }
}

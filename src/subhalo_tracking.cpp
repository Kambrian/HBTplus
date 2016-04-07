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
  if(Nbound>=LastMaxMass) 
  {
	SnapshotIndexOfLastMaxMass=epoch.GetSnapshotIndex();
	LastMaxMass=Nbound;
  }
}
HBTReal Subhalo_t::KineticDistance(const Halo_t &halo, const Snapshot_t &epoch)
{
  HBTReal dx=PeriodicDistance(halo.ComovingAveragePosition, ComovingAveragePosition);
  HBTReal dv=Distance(halo.PhysicalAverageVelocity, PhysicalAverageVelocity);
  HBTReal d=dv+epoch.Hz*epoch.ScaleFactor*dx;
  return (d>0?d:-d);
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
	return (*Subhalos)[i].Nbound>(*Subhalos)[j].Nbound;
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
#ifdef ALLOW_BINARY_SYSTEM
	HBTInt rankoffset=0;
	if(SubGroup.size()>1)
	{
	  if(SubGroup[1]>SubGroup[0]*HBTConfig.BinaryMassRatioLimit)
		rankoffset=1;//binary system, rank start from 1.
	}
#endif
	for(HBTInt i=0;i<SubGroup.size();i++)
#ifdef ALLOW_BINARY_SYSTEM
	  Subhalos[SubGroup[i]].Rank=i+rankoffset;
#else
	  Subhalos[SubGroup[i]].Rank=i;
#endif
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
	Subhalos[subid].HostHaloId=GetLocalHostId(Subhalos[subid].Particles[0].Id, halo_snap, part_snap);
  
  HBTInt nsub=0;
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	if(Subhalos[subid].HostHaloId<0)
	{
	  if(subid>nsub)
		Subhalos[nsub++]=move(Subhalos[subid]);//there should be a default move assignement operator.
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
  vector <SizeRank_t> LocalHostIds(NumSubhalos), GlobalHostIds(NumSubhalos);
  if(thisrank==root)
  {
	#pragma omp parallel for if(NumSubhalos>20)
	for(HBTInt i=0;i<NumSubhalos;i++)
	{
	  LocalHostIds[i].n=Subhalos[i].HostHaloId;//already found previously
	  LocalHostIds[i].rank=thisrank;
	}
  }
  else
  {
	#pragma omp parallel for if(NumSubhalos>20)
	for(HBTInt i=0;i<NumSubhalos;i++)
	{
	  LocalHostIds[i].n=GetLocalHostId(TrackParticleIds[i], halo_snap, part_snap);
	  LocalHostIds[i].rank=thisrank;
	}
  }

  MPI_Allreduce(LocalHostIds.data(), GlobalHostIds.data(), NumSubhalos, MPI_HBTRankPair, MPI_MAXLOC, world.Communicator);

  //scatter free subhaloes from root to everywhere
  //send particles; no scatterw, do it manually
  MPI_Datatype MPI_HBT_Particle;
  Particle_t().create_MPI_type(MPI_HBT_Particle);
  vector <vector<int> > SendSizes(world.size());  
  vector <MPI_Request> Req0,Req1;
  if(thisrank==root)
  {
	vector <vector <MPI_Aint> > SendBuffers(world.size());
	for(HBTInt subid=0;subid<NumSubhalos;subid++)//packing
	{
	  int rank=GlobalHostIds[subid].rank;
	  auto & Particles=Subhalos[subid].Particles;
	  MPI_Aint p;
	  MPI_Address(Particles.data(),&p);
	  SendBuffers[rank].push_back(p);
	  SendSizes[rank].push_back(Particles.size());
	}
	Req0.resize(world.size());
	Req1.resize(world.size());
	for(int rank=0;rank<world.size();rank++)
	{
	  MPI_Isend(SendSizes[rank].data(), SendSizes[rank].size(), MPI_INT, rank, 0, world.Communicator, &Req0[rank]);
	  MPI_Datatype SendType;
	  MPI_Type_create_hindexed(SendSizes[rank].size(), SendSizes[rank].data(), SendBuffers[rank].data(), MPI_HBT_Particle, &SendType);
	  MPI_Type_commit(&SendType);
	  MPI_Isend(MPI_BOTTOM, 1, SendType, rank, 1, world.Communicator, &Req1[rank]);
	  MPI_Type_free(&SendType);
	}
  }
  //receive on every process, including root
	vector <MPI_Aint> ReceiveBuffer;
	vector <int> ReceiveSize;
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
	
	if(thisrank==root)
	{
// 	  vector <MPI_Status> stats(world.size());
	  MPI_Waitall(world.size(), Req0.data(), MPI_STATUS_IGNORE);
	  MPI_Waitall(world.size(), Req1.data(), MPI_STATUS_IGNORE);
	}
  
  //copy other properties
    vector <int> Counts(world.size()), Disps(world.size());
	vector <Subhalo_t> TmpHalos;
	if(world.rank()==root)
	{//reuse GlobalHostIds for sorting
	  for(HBTInt subid=0;subid<Subhalos.size();subid++)
	  {
		Subhalos[subid].HostHaloId=GlobalHostIds[subid].n;
// 		assert(GlobalHostIds[subid].n>=-1);
		GlobalHostIds[subid].n=subid;
	  }
	  stable_sort(GlobalHostIds.begin(), GlobalHostIds.end(), CompareRank);
	  TmpHalos.resize(Subhalos.size());
	  for(HBTInt subid=0;subid<Subhalos.size();subid++)
		TmpHalos[subid]=move(Subhalos[GlobalHostIds[subid].n]);
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
    
  MemberTable.Build(halo_snap.Halos.size(), Subhalos, false);//build without orphans first.
//   MemberTable.AssignRanks(Subhalos); //not needed here
}

void SubhaloSnapshot_t::DecideCentrals(const HaloSnapshot_t &halo_snap)
/* to select central subhalo according to KineticDistance, and move each central to the beginning of each list in MemberTable*/
{
  //initialize the ranks
#ifdef ALLOW_BINARY_SYSTEM
#pragma omp for
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
	Subhalos[subid].Rank=0;
#endif
#pragma omp for
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
  {
	MemberShipTable_t::MemberList_t &List=MemberTable.SubGroups[hostid];
	if(List.size()>1)
	{
#ifdef ALLOW_BINARY_SYSTEM
	  if(Subhalos[List[1]].Nbound>Subhalos[List[0]].Nbound*HBTConfig.BinaryMassRatioLimit)
	  {
		Subhalos[List[0]].Rank=1;//mark as a binary system. do not FeedCentral.
		continue;
	  }
#endif	  
	  int n_major;
	  HBTInt MassLimit=Subhalos[List[0]].Nbound*HBTConfig.MajorProgenitorMassRatio;
	  for(n_major=1;n_major<List.size();n_major++)
		if(Subhalos[List[n_major]].Nbound<MassLimit) break;
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
	if(0==Members.size()) //create a new sub
	{
	  HBTInt subid;
	  #pragma omp critical(AddNewSub) //maybe consider ordered for here..
	  {
		subid=Npro++;
	  }
	  Subhalos[subid].HostHaloId=hostid;
	  copyHBTxyz(Subhalos[subid].ComovingAveragePosition, halo_snap.Halos[hostid].ComovingAveragePosition); 
	  copyHBTxyz(Subhalos[subid].PhysicalAverageVelocity, halo_snap.Halos[hostid].PhysicalAverageVelocity);
	  Subhalos[subid].Particles.swap(halo_snap.Halos[hostid].Particles);
	  
	  Subhalos[subid].SnapshotIndexOfBirth=SnapshotIndex;
	}
	else
	{
#ifdef ALLOW_BINARY_SYSTEM	  
	  if(0==Subhalos[Members[0]].Rank)//only update particles if not a binary system
#endif		
		Subhalos[Members[0]].Particles.swap(halo_snap.Halos[hostid].Particles); //reuse the halo particles
	}
  }
//   #pragma omp single
//   halo_snap.Clear();//to avoid misuse
}
void SubhaloSnapshot_t::PrepareCentrals(HaloSnapshot_t &halo_snap)
{
  #pragma omp parallel
  {
  halo_snap.AverageCoordinates();
  AverageCoordinates();
  DecideCentrals(halo_snap);
  FeedCentrals(halo_snap);
  }
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
void SubhaloSnapshot_t::UpdateTracks(MpiWorker_t &world, const HaloSnapshot_t &halo_snap)
{
  /*renew ranks after unbinding*/
  RegisterNewTracks(world);//performance bottleneck here. no. just poor synchronization.
  #pragma omp parallel
  {
  MemberTable.SortMemberLists(Subhalos);//reorder, so the central might change if necessary
  MemberTable.AssignRanks(Subhalos);
  PurgeMostBoundParticles();
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
  #pragma omp parallel for if(ParallelizeHaloes)
  for(HBTInt i=0;i<Subhalos.size();i++)
  {
	Subhalos[i].CalculateProfileProperties(*this);
	Subhalos[i].CalculateShape();
  }
}

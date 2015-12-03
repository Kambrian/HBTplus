#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"

void Subhalo_t::UpdateTrack(HBTInt snapshot_index)
{
  if(TrackId==SpecialConst::NullTrackId) return;
  
  if(0==Rank) SnapshotIndexOfLastIsolation=snapshot_index;
  if(Nbound>=LastMaxMass) 
  {
	SnapshotIndexOfLastMaxMass=snapshot_index;
	LastMaxMass=Nbound;
  }
}
HBTReal Subhalo_t::KineticDistance(const Halo_t &halo, const ParticleSnapshot_t &snapshot)
{
  HBTReal dx=PeriodicDistance(halo.ComovingPosition, ComovingPosition);
  HBTReal dv=distance(halo.PhysicalVelocity, PhysicalVelocity);
  HBTReal d=dv+snapshot.Hz*snapshot.ScaleFactor*dx;
  return (d>0?d:-d);
}
void MemberShipTable_t::ResizeAllMembers(size_t n)
{
  Mem_AllMembers.resize(n);
  size_t offset=Mem_AllMembers.data()-AllMembers.data();
  if(offset)
  {
	AllMembers.Bind(n, AllMembers.data()+offset);
	for(HBTInt i=0;i<Mem_SubGroups.size();i++)
	  Mem_SubGroups[i].Bind(Mem_SubGroups[i].data()+offset);
  }
}
void MemberShipTable_t::Init(const HBTInt nhalos, const HBTInt nsubhalos, const float alloc_factor)
{
  Mem_SubGroups.clear();
  Mem_SubGroups.resize(nhalos+1);
  SubGroups.Bind(nhalos, Mem_SubGroups.data()+1);

  Mem_AllMembers.clear();
  Mem_AllMembers.reserve(nsubhalos*alloc_factor); //allocate more for seed haloes.
  Mem_AllMembers.resize(nsubhalos);
  AllMembers.Bind(nsubhalos, Mem_AllMembers.data());
}
void MemberShipTable_t::BindMemberLists()
{
  HBTInt offset=0;
  for(HBTInt i=0;i<Mem_SubGroups.size();i++)
  {
	Mem_SubGroups[i].Bind(Mem_SubGroups[i].size(), &(Mem_AllMembers[offset]));
	offset+=Mem_SubGroups[i].size();
	Mem_SubGroups[i].ReBind(0);
  }
}
void MemberShipTable_t::CountMembers(const SubhaloList_t& Subhalos)
{//todo: parallelize this..
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
	SubGroups[Subhalos[subid].HostHaloId].IncrementBind();
}
void MemberShipTable_t::FillMemberLists(const SubhaloList_t& Subhalos)
{//fill with local subhaloid
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
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
	for(HBTInt i=0;i<SubGroup.size();i++)
	  Subhalos[SubGroup[i]].Rank=i;
  }
}
void MemberShipTable_t::CountBirth()
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
void MemberShipTable_t::Build(const HBTInt nhalos, const SubhaloList_t & Subhalos)
{
  #pragma omp single
  {
  Init(nhalos, Subhalos.size());
  CountMembers(Subhalos);
  BindMemberLists();
  FillMemberLists(Subhalos);
  }
  SortMemberLists(Subhalos);
  CountBirth();
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
  
  HBTInt nsub_old=Subhalos.size(), nsub=0;
  for(HBTInt subid=0;subid<nsub;subid++)
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

void FindOtherHosts(mpi::communicator &world, int root, const HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap, vector <Subhalo_t> &Subhalos, vector <Subhalo_t> &LocalSubhalos, MPI_Datatype MPI_Subhalo_Shell_Type)
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
  MPI_Bcast(&NumSubhalos, 1, MPI_HBT_INT, root, world);
  TrackParticleIds.resize(NumSubhalos);
  if(thisrank==root)
  {
	for(HBTInt i=0;i<Subhalos.size();i++)
	  TrackParticleIds[i]=Subhalos[i].Particles[0].Id;
  }
  MPI_Bcast(TrackParticleIds.data(), NumSubhalos, MPI_HBT_INT, root, world);
  
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

  MPI_Allreduce(LocalHostIds.data(), GlobalHostIds.data(), NumSubhalos, MPI_HBTRankPair, MPI_MAXLOC, world);

  //scatter free subhaloes from root to everywhere
  //send particles; no scatterw, do it manually
  MPI_Datatype MPI_HBT_Particle;
  create_MPI_Particle_type(MPI_HBT_Particle);
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
	{//todo: have to use Isend here..... to receive on root later; otherwise root send will deadlock waiting for receive to return.
	  MPI_Isend(SendSizes[rank].data(), SendSizes[rank].size(), MPI_INT, rank, 0, world, &Req0[rank]);
	  MPI_Datatype SendType;
	  MPI_Type_create_hindexed(SendSizes[rank].size(), SendSizes[rank].data(), SendBuffers[rank].data(), MPI_HBT_Particle, &SendType);
	  MPI_Type_commit(&SendType);
	  MPI_Isend(MPI_BOTTOM, 1, SendType, rank, 1, world, &Req1[rank]);
	  MPI_Type_free(&SendType);
	}
  }
  //receive on every process, including root
	vector <MPI_Aint> ReceiveBuffer;
	vector <int> ReceiveSize;
	int NumNewSubs;
	MPI_Status stat;
	MPI_Probe(root, 0, world, &stat);
	MPI_Get_count(&stat, MPI_INT, &NumNewSubs);
	ReceiveSize.resize(NumNewSubs);
	MPI_Recv(ReceiveSize.data(), NumNewSubs, MPI_INT, root, 0, world, &stat);
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
	MPI_Recv(MPI_BOTTOM, 1, ReceiveType, root, 1, world, &stat);
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
	MPI_Scatterv(TmpHalos.data(), Counts.data(), Disps.data(), MPI_Subhalo_Shell_Type, &NewSubhalos[0], NumNewSubs, MPI_Subhalo_Shell_Type, root, world);
}
void SubhaloSnapshot_t::AssignHosts(mpi::communicator &world, HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap)
/* find host haloes for subhaloes, and build MemberTable. Each subhalo is moved to the processor of its host halo, with its HostHaloId set to the local haloid of the host*/
{
  ParallelizeHaloes=halo_snap.NumPartOfLargestHalo<0.1*halo_snap.TotNumberOfParticles;//no dominating objects
  
  //exchange subhaloes according to hosts
  vector <Subhalo_t> LocalSubhalos;
  LocalSubhalos.reserve(Subhalos.size());
  halo_snap.FillParticleHash();
  FindLocalHosts(halo_snap, part_snap, Subhalos, LocalSubhalos);
  for(int rank=0;rank<world.size();rank++)
	FindOtherHosts(world, rank, halo_snap, part_snap, Subhalos, LocalSubhalos, MPI_HBT_SubhaloShell_t);
  Subhalos.swap(LocalSubhalos);
  halo_snap.ClearParticleHash();
    
  MemberTable.Build(halo_snap.Halos.size(), Subhalos);
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
	  HBTInt MassLimit=Subhalos[List[0]].Nbound*HBTConfig.MajorProgenitorMassRatio;
	  for(n_major=1;n_major<List.size();n_major++)
		if(Subhalos[List[n_major]].Nbound<MassLimit) break;
	  if(n_major>1)
	  {
		HBTReal dmin=Subhalos[List[0]].KineticDistance(halo_snap.Halos[hostid], *SnapshotPointer);
		int icenter=0;
		for(int i=1;i<n_major;i++)
		{
		  HBTReal d=Subhalos[List[i]].KineticDistance(halo_snap.Halos[hostid], *SnapshotPointer);
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
	  copyHBTxyz(Subhalos[subid].ComovingPosition, halo_snap.Halos[hostid].ComovingPosition); 
	  copyHBTxyz(Subhalos[subid].PhysicalVelocity, halo_snap.Halos[hostid].PhysicalVelocity);
	  Subhalos[subid].Particles.swap(halo_snap.Halos[hostid].Particles);
	}
	else
	  Subhalos[Members[0]].Particles.swap(halo_snap.Halos[hostid].Particles); //reuse the halo particles
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

void SubhaloSnapshot_t::RegisterNewTracks(mpi::communicator &world)
/*assign trackId to new bound ones, remove unbound ones, and record membership*/
{
  HBTInt NTot=Subhalos.size();
  HBTInt Nsub=NTot-MemberTable.NBirth;
  MemberTable.ResizeAllMembers(NTot);
  for(HBTInt i=Nsub;i<NTot;i++)
  {
	if(Subhalos[i].Nbound>1)
	{
	  if(i!=Nsub)
		Subhalos[Nsub]=move(Subhalos[i]);
	  MemberTable.AllMembers[Nsub]=Nsub; //the MemberTable stores local subid, not TrackId.
	  MemberTable.SubGroups[Subhalos[Nsub].HostHaloId].Bind(1, &MemberTable.AllMembers[Nsub]);
	  Nsub++;
	}
  }
  MemberTable.NFake=NTot-Nsub;
  MemberTable.NBirth-=MemberTable.NFake;
  Subhalos.resize(Nsub);
//   MemberTable.ResizeAllMembers(Nsub); //not necessary
  
  //now assign a global TrackId
  HBTInt TrackIdOffset, NBirth=MemberTable.NBirth, Nsub_last=Nsub-NBirth, GlobalNumberOfSubs;
  MPI_Allreduce(&Nsub_last, &GlobalNumberOfSubs, 1, MPI_HBT_INT, MPI_SUM, world); 
  MPI_Scan(&NBirth, &TrackIdOffset, 1, MPI_HBT_INT, MPI_SUM, world); 
  TrackIdOffset=TrackIdOffset+GlobalNumberOfSubs-NBirth;
  for(HBTInt i=Nsub_last;i<Nsub;i++)
	Subhalos[i].TrackId=TrackIdOffset++;
}
void SubhaloSnapshot_t::UpdateTracks(mpi::communicator &world, const HaloSnapshot_t &halo_snap)
{
  /*renew ranks after unbinding*/
  RegisterNewTracks(world);
  #pragma omp parallel
  {
  MemberTable.SortMemberLists(Subhalos);//reorder, so the central might change if necessary
  MemberTable.AssignRanks(Subhalos);
  #pragma omp for
  for(HBTInt i=0;i<Subhalos.size();i++)
  {
	Subhalos[i].UpdateTrack(SnapshotIndex);
	HBTInt HostId=Subhalos[i].HostHaloId;
	if(HostId<0)
	  Subhalos[i].HostHaloId=-1;
	else
	  Subhalos[i].HostHaloId=halo_snap.Halos[HostId].HaloId;//restore global haloid
  }
  }
}

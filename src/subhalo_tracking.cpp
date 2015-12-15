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
{
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

void SubhaloSnapshot_t::AssignHosts(const HaloSnapshot_t &halo_snap)
{
  static vector<HBTInt> ParticleToHost;//to make it shared
#pragma omp single
{
   ParticleToHost.assign(SnapshotPointer->GetNumberOfParticles(), SpecialConst::NullHaloId);
   ParallelizeHaloes=halo_snap.NumPartOfLargestHalo<0.1*halo_snap.TotNumberOfParticles;//no dominating objects
}
#pragma omp for
  for(HBTInt haloid=0;haloid<halo_snap.Halos.size();haloid++)
  {
	const Halo_t::ParticleList_t & Particles=halo_snap.Halos[haloid].Particles;
	for(HBTInt i=0;i<Particles.size();i++)
	  ParticleToHost[Particles[i]]=haloid;
  }
#pragma omp for 
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	//rely on most-bound particle
	Subhalos[subid].HostHaloId=ParticleToHost[Subhalos[subid].Particles[0]];
	//alternatives: CoreTrack; Split;
  }
#pragma omp single
  vector <HBTInt>().swap(ParticleToHost);//free the memory.
#pragma omp single nowait
  if(HBTConfig.TrimNonHostParticles)
  {
	cout<<"Error: TrimNonHostParticles not implemented yet...\n";
	exit(1);
  }
  MemberTable.Build(halo_snap.Halos.size(), Subhalos);
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
  //Subhalos.reserve(Snapshot->GetNumberOfParticles()*0.1);//reserve enough	branches.......
  Subhalos.resize(Npro+MemberTable.NBirth);
  }
  #pragma omp for
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
  { 
	MemberShipTable_t::MemberList_t & Members=MemberTable.SubGroups[hostid];
	if(0==Members.size()) //create a new sub
	{
	  HBTInt subid;
	  #pragma omp critical(AddNewSub) //maybe consider ordered for here, to ensure group order?
	  {
		subid=Npro++;
	  }
	  Subhalos[subid].HostHaloId=hostid;
	  copyHBTxyz(Subhalos[subid].ComovingPosition, halo_snap.Halos[hostid].ComovingPosition); 
	  copyHBTxyz(Subhalos[subid].PhysicalVelocity, halo_snap.Halos[hostid].PhysicalVelocity);
	  Subhalos[subid].Particles.swap(halo_snap.Halos[hostid].Particles);
	}
	else
	{
#ifdef ALLOW_BINARY_SYSTEM	  
	  if(0==Subhalos[Members[0]].Rank)//only update particles if not a binary system
#endif		
		Subhalos[Members[0]].Particles.swap(halo_snap.Halos[hostid].Particles); //reuse the halo particles
	}
  }
  #pragma omp single
  halo_snap.Clear();//to avoid misuse
}
void SubhaloSnapshot_t::PrepareCentrals(HaloSnapshot_t &halo_snap)
{
  halo_snap.AverageCoordinates();
  AverageCoordinates();
  DecideCentrals(halo_snap);
  FeedCentrals(halo_snap);
}

void SubhaloSnapshot_t::RegisterNewTracks()
/*assign trackId to new bound ones, remove unbound ones, and record membership*/
{
  HBTInt NTot=Subhalos.size();
  HBTInt TrackId=NTot-MemberTable.NBirth, NFake=0;
  MemberTable.ResizeAllMembers(NTot);
  for(HBTInt i=TrackId;i<NTot;i++)
  {
	if(Subhalos[i].Nbound>1)
	{
	  if(i!=TrackId)
		Subhalos[i].MoveTo(Subhalos[TrackId]);
	  Subhalos[TrackId].TrackId=TrackId;
	  MemberTable.AllMembers[TrackId]=TrackId; //this trackId is also the subhalo index
// 	  assert(MemberTable.SubGroups[Subhalos[TrackId].HostHaloId].size()==0);
	  MemberTable.SubGroups[Subhalos[TrackId].HostHaloId].Bind(1, &MemberTable.AllMembers[TrackId]);
	  TrackId++;
	}
  }
  MemberTable.NFake=NTot-TrackId;
  MemberTable.NBirth-=MemberTable.NFake;
  Subhalos.resize(TrackId);
}
void SubhaloSnapshot_t::UpdateTracks()
{
  /*renew ranks after unbinding*/
#pragma omp single
  RegisterNewTracks();
  MemberTable.SortMemberLists(Subhalos);//reorder
  MemberTable.AssignRanks(Subhalos);
#pragma omp for
  for(HBTInt i=0;i<Subhalos.size();i++)
	Subhalos[i].UpdateTrack(SnapshotIndex);
}

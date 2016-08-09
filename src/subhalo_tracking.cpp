#include <iostream>
#include <new>
#include <algorithm>
#include <unordered_set>
#include <omp.h>

#include "datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"

void Subhalo_t::UpdateTrack(const ParticleSnapshot_t &part_snap)
{
  if(TrackId==SpecialConst::NullTrackId) return;
  
  if(0==Rank) SnapshotIndexOfLastIsolation=part_snap.GetSnapshotIndex();
  if(Mbound>=LastMaxMass) 
  {
	SnapshotIndexOfLastMaxMass=part_snap.GetSnapshotIndex();
	LastMaxMass=Mbound;
  }
}
HBTReal Subhalo_t::KineticDistance(const Halo_t &halo, const ParticleSnapshot_t &snapshot)
{
  HBTxyz dv;
  snapshot.RelativeVelocity(ComovingAveragePosition, PhysicalAverageVelocity, halo.ComovingAveragePosition, halo.PhysicalAverageVelocity, dv);
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
{
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
  
  if(a.HostHaloId==b.HostHaloId) return (a.Mbound>b.Mbound);
  
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
	auto &Particles=Subhalos[subid].Particles;
	if(Particles.size())
	  Subhalos[subid].HostHaloId=ParticleToHost[Subhalos[subid].Particles[0]];
	else
	  Subhalos[subid].HostHaloId=SpecialConst::NullHaloId;
  }
#pragma omp single
  vector <HBTInt>().swap(ParticleToHost);//free the memory.
/*#pragma omp single nowait
  if(HBTConfig.TrimNonHostParticles)
  {
	cout<<"Error: TrimNonHostParticles not implemented yet...\n";
	exit(1);
  } */
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
	  HBTReal MassLimit=Subhalos[List[0]].Mbound*HBTConfig.MajorProgenitorMassRatio;
	  for(n_major=1;n_major<List.size();n_major++)
		if(Subhalos[List[n_major]].Mbound<MassLimit) break;
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
	auto &Host=halo_snap.Halos[hostid];
	MemberShipTable_t::MemberList_t & Members=MemberTable.SubGroups[hostid];
	if(0==Members.size()) //create a new sub
	{
	  HBTInt subid;
	  #pragma omp critical(AddNewSub) //maybe consider ordered for here, to ensure group order?
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
	  {
		central.Particles.swap(Host.Particles); //reuse the halo particles
		central.Nbound=central.Particles.size();
		auto &mostbndid=Host.Particles[0]; 
		for(auto & p: central.Particles)
		if(p==mostbndid)//swap previous mostbound particle to the beginning
		{
		  swap(p, central.Particles[0]);
		  break;
		}
	  }
	}
  }
  #pragma omp single
  halo_snap.Clear();//to avoid misuse
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
  NestSubhalos();
#ifndef INCLUSIVE_MASS
  MaskSubhalos();
#endif
}

void SubhaloSnapshot_t::RegisterNewTracks()
/*assign trackId to new bound ones, remove unbound ones, and rebuild membership*/
{
  HBTInt NumSub=Subhalos.size(), NumTrackOld=NumSub-MemberTable.NBirth;
  HBTInt TrackId=NumTrackOld;
  for(HBTInt i=TrackId;i<NumSub;i++)
  {
	if(Subhalos[i].Nbound>1)
	{
	  if(i!=TrackId)
		Subhalos[TrackId]=move(Subhalos[i]);
	  Subhalos[TrackId].TrackId=TrackId;
	  TrackId++;
	}
  }
  Subhalos.resize(TrackId);//erase unbound ones
  MemberTable.Build(MemberTable.SubGroups.size(), Subhalos, true);//rebuild membership with new subs and also include orphans this time.
  MemberTable.NFake=NumSub-TrackId;
  MemberTable.NBirth=TrackId-NumTrackOld;
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
		  if(ExclusionList.find(p)==ExclusionList.end())
		  {
			if(&p!=&subhalo.Particles[0])
			{
			  copyHBTxyz(subhalo.ComovingMostBoundPosition, SnapshotPointer->GetComovingPosition(p));
			  copyHBTxyz(subhalo.PhysicalMostBoundVelocity, SnapshotPointer->GetPhysicalVelocity(p));
			  swap(subhalo.Particles[0], p);
			  break;
			}
		  }
		}
		ExclusionList.insert(subhalo.Particles.begin(), subhalo.Particles.end());//alternative: only filter most-bounds
	  }
	}
  }
}

void SubhaloSnapshot_t::NestSubhalos()
{
  LevelUpDetachedSubhalos();
  //collect detached(head) subhalos
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
  vector <bool> IsHeadSub(Subhalos.size());
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
      auto insert_status=ExclusionList.insert(*it);
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

void SubhaloSnapshot_t::UpdateTracks()
{
  /*renew ranks after unbinding*/
  RegisterNewTracks();
#pragma omp parallel
  {
  MemberTable.SortMemberLists(Subhalos);//reorder
  ExtendCentralNest();
  MemberTable.AssignRanks(Subhalos);
#ifdef INCLUSIVE_MASS
  PurgeMostBoundParticles();
#endif
#pragma omp for
  for(HBTInt i=0;i<Subhalos.size();i++)
	Subhalos[i].UpdateTrack(*SnapshotPointer);
  }
#pragma omp parallel for if(ParallelizeHaloes)
  for(HBTInt i=0;i<Subhalos.size();i++)
  {
	Subhalos[i].CalculateProfileProperties(*SnapshotPointer);
	Subhalos[i].CalculateShape(*SnapshotPointer);
  }
}
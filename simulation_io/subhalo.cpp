#include <iostream>
#include <new>
#include <algorithm>

#include "../datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"
#include "../gravity/tree.h"

HBTReal SubHalo_t::KineticDistance(const Halo_t &halo, const Snapshot_t &snapshot)
{
  HBTReal dx=PeriodicDistance(halo.ComovingPosition, ComovingPosition);
  HBTReal dv=distance(halo.PhysicalVelocity, PhysicalVelocity);
  HBTReal d=dv+snapshot.Header.Hz*snapshot.Header.time*dx;
  return (d>0?d:-d);
}
void SubHalo_t::unbind(const Snapshot_t &snapshot)
{
  OctTree_t tree;
  tree.Build(Particles.size(), Particles.data(), &(snapshot.GetComovingPosition(0)));
  /*FIXME: finish this*/
}

void MemberShipTable_t::Init(const HBTInt nhalos, const HBTInt nsubhalos, const float alloc_factor)
{
  RawLists.clear();
  RawLists.resize(nhalos+1);
  Lists.Bind(nhalos, RawLists.data()+1);

  AllMembers.clear();
  AllMembers.reserve(nsubhalos*alloc_factor); //allocate more for seed haloes.
  AllMembers.resize(nsubhalos);
}
void MemberShipTable_t::BindMemberLists()
{
  HBTInt offset=0;
  for(HBTInt i=0;i<RawLists.size();i++)
  {
	RawLists[i].Bind(RawLists[i].size(), &(AllMembers[offset]));
	offset+=RawLists[i].size();
	RawLists[i].Resize(0);
  }
}
void MemberShipTable_t::CountMembers(const SubHaloList_t& SubHalos)
{
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
	Lists[SubHalos[subid].HostHaloId].IncrementSize();
}
void MemberShipTable_t::FillMemberLists(const SubHaloList_t& SubHalos)
{
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
	Lists[SubHalos[subid].HostHaloId].push_back(subid);
}
struct CompareMass_t
{
  const SubHaloList_t * SubHalos;
  CompareMass_t(const SubHaloList_t & subhalos)
  {
	SubHalos=&subhalos;
  }
  bool operator () (const HBTInt & i, const HBTInt & j)
  {
	return (*SubHalos)[i].Nbound>(*SubHalos)[j].Nbound;
  }
};
void MemberShipTable_t::SortMemberLists(const SubHaloList_t & SubHalos)
{
  CompareMass_t compare_mass(SubHalos);
  for(HBTInt i=0;i<RawLists.size();i++)
	std::sort(RawLists[i].data(), RawLists[i].data()+RawLists[i].size(), compare_mass);
}
void MemberShipTable_t::Build(const HBTInt nhalos, const SubHaloList_t & SubHalos)
{
  Init(nhalos, SubHalos.size());
  CountMembers(SubHalos);
  BindMemberLists();
  FillMemberLists(SubHalos);
  SortMemberLists(SubHalos);
//   std::sort(AllMembers.begin(), AllMembers.end(), CompareHostAndMass);
}
/*
inline bool SubHaloSnapshot_t::CompareHostAndMass(const HBTInt& subid_a, const HBTInt& subid_b)
{//ascending in host id, descending in mass inside each host, and put NullHaloId to the beginning.
  SubHalo_t a=SubHalos[subid_a], b=SubHalos[subid_b];
  
  if(a.HostHaloId==b.HostHaloId) return (a.Nbound>b.Nbound);
  
  return (a.HostHaloId<b.HostHaloId); //(a.HostHaloId!=SpecialConst::NullHaloId)&&
}*/
void SubHaloSnapshot_t::AssignHost(const HaloSnapshot_t &halo_snap)
{
  vector<HBTInt> ParticleToHost(SnapshotPointer->GetNumberOfParticles(), SpecialConst::NullHaloId);
  for(int haloid=0;haloid<halo_snap.Halos.size();haloid++)
  {
	for(int i=0;i<halo_snap.Halos[haloid].Particles.size();i++)
	  ParticleToHost[i]=haloid;
  }
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
  {
	//rely on track_particle
	SubHalos[subid].HostHaloId=ParticleToHost[SubHalos[subid].TrackParticleId];
	//alternatives: CoreTrack; Split;
  }
  //alternative: trim particles outside fof
    
  MemberTable.Build(halo_snap.Halos.size(), SubHalos);
}
#define CORE_SIZE_MIN 20
#define CORE_SIZE_FRAC 0.25
void SubHaloSnapshot_t::AverageCoordinates()
{
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
  {
	int coresize=SubHalos[subid].Nbound*CORE_SIZE_FRAC;
	if(coresize<CORE_SIZE_MIN) coresize=CORE_SIZE_MIN;
	if(coresize>SubHalos[subid].Nbound) coresize=SubHalos[subid].Nbound;
	
	SnapshotPointer->AveragePosition(SubHalos[subid].ComovingPosition, SubHalos[subid].Particles.data(), coresize);
	SnapshotPointer->AverageVelocity(SubHalos[subid].PhysicalVelocity, SubHalos[subid].Particles.data(), coresize);
  }
}

#define MAJOR_PROGENITOR_MASS_RATIO 0.67
void SubHaloSnapshot_t::DecideCentrals(const HaloSnapshot_t &halo_snap)
/* to select central subhalo according to KineticDistance, and move each central to the beginning of each list in MemberTable*/
{
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
  {
	MemberShipTable_t::MemberList_t &List=MemberTable.Lists[hostid];
	if(List.size()>1)
	{
	  int n_major;
	  HBTInt MassLimit=SubHalos[List[0]].Nbound*MAJOR_PROGENITOR_MASS_RATIO;
	  for(n_major=1;n_major<List.size();n_major++)
		if(SubHalos[List[n_major]].Nbound<MassLimit) break;
		if(n_major>1)
		{
		  HBTReal dmin=SubHalos[List[0]].KineticDistance(halo_snap.Halos[hostid], *SnapshotPointer);
		  int icenter=0;
		  for(int i=1;i<n_major;i++)
		  {
			HBTReal d=SubHalos[List[i]].KineticDistance(halo_snap.Halos[hostid], *SnapshotPointer);
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

void SubHaloSnapshot_t::FeedCentrals(HaloSnapshot_t& halo_snap)
{
  HBTInt Npro=SubHalos.size(),NBirth=0;
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
	if(MemberTable.Lists[hostid].size()==0)  NBirth++;
  SubHalos.resize(Npro+NBirth);
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
  {
	if(0==MemberTable.Lists[hostid].size()) //create a new sub
	{//assign trackId only at save time.
		SubHalos[Npro].HostHaloId=hostid;
		SubHalos[Npro].SnapshotIndexOfBirth=SnapshotIndex;
		MemberTable.AllMembers.push_back(Npro);
		MemberTable.Lists[hostid].Bind(1, &MemberTable.AllMembers.back());
		Npro++;
	}
	SubHalos[MemberTable.Lists[hostid][0]].Particles.swap(halo_snap.Halos[hostid].Particles); //reuse the halo particles; now halo_snap should
  }
}


void SubHaloSnapshot_t::RefineParticles()
{//it's more expensive to build an exclusive list. so do inclusive here. 
  //TODO: ensure the inclusive unbinding is stable (contaminating particles from big subhaloes may hurdle the unbinding
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
  {
	SubHalos[subid].unbind(*SnapshotPointer);
  }
}

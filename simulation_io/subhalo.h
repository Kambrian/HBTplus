#ifndef SUBHALO_HEADER_INCLUDED
#define SUBHALO_HEADER_INCLUDED

#include <iostream>
#include <new>
#include <vector>

#include "../datatypes.h"
#include "snapshot_number.h"
#include "halo.h"

class TrackParticle_t
{
public:
  HBTInt TrackId;
  HBTInt SnapshotIndexOfLastIsolation; //the last snapshot when it was a central, only considering past snapshots.
  HBTInt SnapshotIndexOfLastMaxMass; //the snapshot when it has the maximum subhalo mass, only considering past snapshots.
//   HBTInt SnapshotIndexOfBirth;//when the subhalo first becomes resolved
//   HBTInt SnapshotIndexOfDeath;//when the subhalo first becomes un-resolved.
  HBTInt LastMaxMass;
  TrackParticle_t()
  {
	TrackId=SpecialConst::NullTrackId;
	SnapshotIndexOfLastIsolation=SpecialConst::NullSnapshotId;
	SnapshotIndexOfLastMaxMass=SpecialConst::NullSnapshotId;
// 	SnapshotIndexOfBirth=SpecialConst::NullSnapshotId;
// 	SnapshotIndexOfDeath=SpecialConst::NullSnapshotId;
	LastMaxMass=0;
  }
};

class SubHalo_t: public Halo_t, public TrackParticle_t
{
public:
  HBTInt Nbound;
  HBTInt HostHaloId;
  HBTInt Rank;
//   HBTReal RmaxComoving;
//   HBTReal VmaxPhysical;
//   HBTReal RPoissonComoving;
  SubHalo_t(): Halo_t(), TrackParticle_t(), Nbound(0), Rank(0)
  {
  }
  void Unbind(const Snapshot_t &part_snap);
  HBTReal KineticDistance(const Halo_t & halo, const Snapshot_t & partsnap);
  void UpdateTrack(HBTInt snapshot_index);
  bool IsCentral()
  {
	return 0==Rank;
  }
};

typedef vector <SubHalo_t> SubHaloList_t;

class MemberShipTable_t
/* list the subhaloes inside each host, rather than ordering the subhaloes 
 * 
 * the principle is to not move the objects, but construct a table of them, since moving objects will change their id (or index at least), introducing the trouble to re-index them and update the indexes in any existence references.
 */
{
public:
  typedef ShallowList_t<HBTInt> MemberList_t;  //list of members in a group
private:
  void Init(const HBTInt nhalos, const HBTInt nsubhalos, const float alloc_factor=1.2);
  void BindMemberLists();
  void FillMemberLists(const SubHaloList_t & SubHalos);
  void CountMembers(const SubHaloList_t & SubHalos);
  void SortSatellites(const SubHaloList_t & SubHalos);
  void CountBirth();
  vector <MemberList_t> RawSubGroups; //list of subhaloes inside each host halo; contain one more group than halo catalogue, to hold field subhaloes.
public:
  vector <HBTInt> AllMembers; //the storage for all the MemberList_t
  ShallowList_t <MemberList_t> SubGroups; //list of subhaloes inside each host halo. offset from RawSubGroups to allow hostid=-1.
  HBTInt NBirth; //newly born halos.
  HBTInt NBirthFake; //Fake (unbound) halos out of NBirth
  
  MemberShipTable_t(): RawSubGroups(), AllMembers(), SubGroups(), NBirth(0), NBirthFake(0)
  {
  }
  void Build(const HBTInt nhalos, const SubHaloList_t & SubHalos);
  void SortMemberLists(const SubHaloList_t & SubHalos);
  void AssignRanks(SubHaloList_t &SubHalos);
};
class SubHaloSnapshot_t: public SnapshotNumber_t
{ 
private:
  void RegisterNewTracks();
  void DecideCentrals(const HaloSnapshot_t &halo_snap);
  void FeedCentrals(HaloSnapshot_t &halo_snap);
public:
  Snapshot_t * SnapshotPointer;
  SubHaloList_t SubHalos;
  MemberShipTable_t MemberTable;
  SubHaloSnapshot_t(): SnapshotNumber_t(), SubHalos(), MemberTable(), SnapshotPointer(nullptr)
  {
  }
  void Load(int snapshot_index);
  void Save();
  void Clear()
  {
	//TODO
	cout<<"Clean() not implemented yet\n";
  }
  void ParticleIdToIndex(Snapshot_t & snapshot)
  {//also bind to snapshot
	SnapshotPointer=&snapshot;
	for(HBTInt subid=0;subid<SubHalos.size();subid++)
	for(HBTInt pid=0;pid<SubHalos[pid].Particles.size();pid++)
	  SubHalos[subid].Particles[pid]=snapshot.GetParticleIndex(SubHalos[subid].Particles[pid]);
  }
  void ParticleIndexToId()
  {
	Snapshot_t &snapshot=*SnapshotPointer;
	for(HBTInt subid=0;subid<SubHalos.size();subid++)
	for(HBTInt pid=0;pid<SubHalos[pid].Particles.size();pid++)
	  SubHalos[subid].Particles[pid]=snapshot.GetParticleId(SubHalos[subid].Particles[pid]);
	SnapshotPointer=nullptr;
  }
  void AverageCoordinates();
  void AssignHosts(const HaloSnapshot_t &halo_snap);
  void PrepareCentrals(HaloSnapshot_t &halo_snap);
  void RefineParticles();
  void UpdateRanks();
  void UpdateTracks();
};

#endif
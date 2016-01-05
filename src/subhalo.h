#ifndef SUBHALO_HEADER_INCLUDED
#define SUBHALO_HEADER_INCLUDED

#include <iostream>
#include <new>
#include <vector>
#include "hdf5.h"
#include "hdf5_hl.h"	
// #include "H5Cpp.h"
#ifdef HBT_REAL8
#define H5T_HBTReal H5T_NATIVE_DOUBLE
#else
#define H5T_HBTReal H5T_NATIVE_FLOAT
#endif
#ifdef HBT_INT8
#define H5T_HBTInt H5T_NATIVE_LONG
#else 
#define H5T_HBTInt H5T_NATIVE_INT
#endif

#include "datatypes.h"
#include "snapshot_number.h"
#include "halo.h"

class Subhalo_t
{
public:
  typedef vector <HBTInt> ParticleList_t;
  HBTInt TrackId;
  HBTInt Nbound;
  HBTInt HostHaloId;
  HBTInt Rank;
  HBTInt LastMaxMass;
  HBTInt SnapshotIndexOfLastMaxMass; //the snapshot when it has the maximum subhalo mass, only considering past snapshots.
  HBTInt SnapshotIndexOfLastIsolation; //the last snapshot when it was a central, only considering past snapshots.
  
  HBTInt SnapshotIndexOfBirth;//when the subhalo first becomes resolved
  HBTInt SnapshotIndexOfDeath;//when the subhalo first becomes un-resolved; only set if currentsnapshot>=SnapshotIndexOfDeath.
  
  //profile properties
  HBTReal RmaxComoving;
  HBTReal VmaxPhysical;
  HBTReal LastMaxVmaxPhysical;
  HBTInt SnapshotIndexOfLastMaxVmax; //the snapshot when it has the maximum Vmax, only considering past snapshots.
  
  HBTReal R2SigmaComoving; //95.5% containment radius
  HBTReal RHalfComoving;
  
  HBTReal R200CritComoving;
  HBTReal R200MeanComoving;
  HBTReal RVirComoving;
  HBTReal M200Crit;
  HBTReal M200Mean;
  HBTReal MVir;
  
  //kinetic properties
  HBTReal SpecificSelfPotentialEnergy;
  HBTReal SpecificSelfKineticEnergy;
  HBTReal SpecificAngularMomentum;
  HBTReal SpinPeebles;
  HBTReal SpinBullock;
  
  //shapes
  HBTxyz InertialEigenVector[3];
  HBTxyz InertialEigenVectorWeighted[3];
  
  HBTxyz ComovingAveragePosition;
  HBTxyz PhysicalAverageVelocity;
  HBTxyz ComovingMostBoundPosition;
  HBTxyz PhysicalMostBoundVelocity;
  
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  
  ParticleList_t Particles;
  
  Subhalo_t(): Nbound(0), Rank(0)
  {
	TrackId=SpecialConst::NullTrackId;
	SnapshotIndexOfLastIsolation=SpecialConst::NullSnapshotId;
	SnapshotIndexOfLastMaxMass=SpecialConst::NullSnapshotId;
	LastMaxMass=0;
	
	SnapshotIndexOfBirth=-1;
	SnapshotIndexOfDeath=-1;
  }
  /*void MoveTo(Subhalo_t & dest)
  {//override dest with this, leaving this unspecified.
	dest.TrackId=TrackId;
	dest.Nbound=Nbound;
	dest.HostHaloId=HostHaloId;
	dest.Rank=Rank;
	dest.LastMaxMass=LastMaxMass;
	dest.SnapshotIndexOfLastMaxMass=SnapshotIndexOfLastMaxMass;
	dest.SnapshotIndexOfLastIsolation=SnapshotIndexOfLastIsolation;
	copyHBTxyz(dest.ComovingPosition, ComovingPosition);
	copyHBTxyz(dest.PhysicalVelocity, PhysicalVelocity);
	dest.Particles.swap(Particles);
  }*/
  void Unbind(const ParticleSnapshot_t &part_snap);
  HBTReal KineticDistance(const Halo_t & halo, const ParticleSnapshot_t & partsnap);
  void UpdateTrack(const ParticleSnapshot_t &part_snap);
  bool IsCentral()
  {
	return 0==Rank;
  }
  void CalculateProfileProperties(const ParticleSnapshot_t &part_snap);
};

typedef vector <Subhalo_t> SubhaloList_t;

class MemberShipTable_t
/* list the subhaloes inside each host, rather than ordering the subhaloes 
 * 
 * the principle is to not move the objects, but construct a table of them, since moving objects will change their id (or index at least), introducing the trouble to re-index them and update the indexes in any existence references.
 */
{
public:
  typedef VectorView_t <HBTInt> MemberList_t;  //list of members in a group
private:
  void BindMemberLists();
  void FillMemberLists(const SubhaloList_t & Subhalos);
  void CountMembers(const SubhaloList_t & Subhalos);
  void SortSatellites(const SubhaloList_t & Subhalos);
  void CountBirth();
  /*avoid operating on the Mem_* below; use the public VectorViews whenever possible; only operate the Mem_* variables when adjusting memory*/
  vector <MemberList_t> Mem_SubGroups; //list of subhaloes inside each host halo, with the storage of each subgroup mapped to a location in Mem_AllMembers 
  vector <HBTInt> Mem_AllMembers; //the storage for all the MemberList_t
public:
  VectorView_t <HBTInt> AllMembers; //the complete list of all the subhaloes in SubGroups.
  VectorView_t <MemberList_t> SubGroups; //list of subhaloes inside each host halo. contain one more group than halo catalogue, to hold field subhaloes. It is properly offseted so that SubGroup[hostid=-1] gives field subhaloes, and hostid>=0 for the normal groups.
  HBTInt NBirth; //newly born halos, excluding fake halos
  HBTInt NFake; //Fake (unbound) halos with no progenitors
  
  MemberShipTable_t(): Mem_SubGroups(), Mem_AllMembers(), AllMembers(), SubGroups(), NBirth(0), NFake(0)
  {
  }
  HBTInt GetNumberOfFieldSubs()
  {
	return SubGroups[-1].size();
  }
  void Init(const HBTInt nhalos, const HBTInt nsubhalos, const float alloc_factor=1.2);
  void ResizeAllMembers(size_t n);
  void Build(const HBTInt nhalos, const SubhaloList_t & Subhalos);
  void SortMemberLists(const SubhaloList_t & Subhalos);
  void AssignRanks(SubhaloList_t &Subhalos);
  void SubIdToTrackId(const SubhaloList_t &Subhalos);
  void TrackIdToSubId(SubhaloList_t &Subhalos);
};
class SubhaloSnapshot_t: public Snapshot_t
{ 
private:
  void RegisterNewTracks();
  bool ParallelizeHaloes;
  hid_t H5T_SubhaloInMem, H5T_SubhaloInDisk;
  
  void DecideCentrals(const HaloSnapshot_t &halo_snap);
  void FeedCentrals(HaloSnapshot_t &halo_snap);
  void BuildHDFDataType();
public:
  const ParticleSnapshot_t * SnapshotPointer;
  SubhaloList_t Subhalos;
  MemberShipTable_t MemberTable;
  SubhaloSnapshot_t(): Snapshot_t(), Subhalos(), MemberTable(), SnapshotPointer(nullptr),  ParallelizeHaloes(true)
  {
	BuildHDFDataType();
  }
  void GetSubFileName(string &filename);
  void GetSrcFileName(string &filename);
  void Load(int snapshot_index, bool load_src=false);
  void Save();
  ~SubhaloSnapshot_t()
  {
	H5Tclose(H5T_SubhaloInDisk);
	H5Tclose(H5T_SubhaloInMem);
  }
  void Clear()
  {
	//TODO
	cout<<"Clean() not implemented yet\n";
  }
  void ParticleIdToIndex(const ParticleSnapshot_t & snapshot);
  void ParticleIndexToId();
  void AverageCoordinates();
  void AssignHosts(const HaloSnapshot_t &halo_snap);
  void PrepareCentrals(HaloSnapshot_t &halo_snap);
  void RefineParticles();
  void UpdateTracks();
  HBTInt size() const
  {
	return Subhalos.size();
  }
  HBTInt GetMemberId(const HBTInt index)
  {
	return Subhalos[index].TrackId;
  }
  const HBTxyz & GetComovingPosition(const HBTInt index) const
  {
	return Subhalos[index].ComovingPosition;
  }
  const HBTxyz & GetPhysicalVelocity(const HBTInt index) const
  {
	return Subhalos[index].PhysicalVelocity;
  }
  HBTReal GetMass(const HBTInt index) const
  {
	return Subhalos[index].Particles.size();
  }
};
inline HBTInt GetCoreSize(HBTInt nbound)
/* get the size of the core that determines the position of the subhalo.
 * coresize controlled by SubCoreSizeFactor and SubCoreSizeMin.
 * if you do not want a cored center, then
 * set SubCoreSizeFactor=0 and SubCoreSizeMin=1 to use most-bound particle;
 * set SubCoreSizeFactor=1 to use all the particles*/
{
  int coresize=nbound*HBTConfig.SubCoreSizeFactor;
  if(coresize<HBTConfig.SubCoreSizeMin) coresize=HBTConfig.SubCoreSizeMin;
  if(coresize>nbound) coresize=nbound;
  return coresize;
}
#endif
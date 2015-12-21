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
  typedef vector <Particle_t> ParticleList_t;
  HBTInt TrackId;
  HBTInt Nbound;
  HBTInt HostHaloId;
  HBTInt Rank;
  HBTInt LastMaxMass;
  HBTInt SnapshotIndexOfLastMaxMass; //the snapshot when it has the maximum subhalo mass, only considering past snapshots.
  HBTInt SnapshotIndexOfLastIsolation; //the last snapshot when it was a central, only considering past snapshots.
//   HBTInt SnapshotIndexOfBirth;//when the subhalo first becomes resolved
//   HBTInt SnapshotIndexOfDeath;//when the subhalo first becomes un-resolved.
//   HBTReal RmaxComoving;
//   HBTReal VmaxPhysical;
//   HBTReal RPoissonComoving;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  ParticleList_t Particles;
  
  Subhalo_t(): Nbound(0), Rank(0)
  {
	TrackId=SpecialConst::NullTrackId;
	SnapshotIndexOfLastIsolation=SpecialConst::NullSnapshotId;
	SnapshotIndexOfLastMaxMass=SpecialConst::NullSnapshotId;
	LastMaxMass=0;
  }
 /* deprecated. use move assignement instead.
  * void MoveTo(Subhalo_t & dest, bool MoveParticle=true)
  {//override dest with this, leaving this unspecified if MoveParticle=true.
	dest.TrackId=TrackId;
	dest.Nbound=Nbound;
	dest.HostHaloId=HostHaloId;
	dest.Rank=Rank;
	dest.LastMaxMass=LastMaxMass;
	dest.SnapshotIndexOfLastMaxMass=SnapshotIndexOfLastMaxMass;
	dest.SnapshotIndexOfLastIsolation=SnapshotIndexOfLastIsolation;
	copyHBTxyz(dest.ComovingPosition, ComovingPosition);
	copyHBTxyz(dest.PhysicalVelocity, PhysicalVelocity);
	if(MoveParticle)
	  dest.Particles.swap(Particles);
  }*/
  void Unbind(const Snapshot_t &epoch);
  HBTReal KineticDistance(const Halo_t & halo, const Snapshot_t & partsnap);
  void UpdateTrack(HBTInt snapshot_index);
  bool IsCentral()
  {
	return 0==Rank;
  }
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
  VectorView_t <HBTInt> AllMembers; //the complete list of all the subhaloes in SubGroups. (contains local subhaloid, i.e., the index of subhaloes in the local SubhaloSnapshot_t)
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
  bool ParallelizeHaloes;
  hid_t H5T_SubhaloInMem, H5T_SubhaloInDisk;
  MPI_Datatype MPI_HBT_SubhaloShell_t;//MPI datatype ignoring the particle list
  
  void RegisterNewTracks(mpi::communicator &world);
  void DecideCentrals(const HaloSnapshot_t &halo_snap);
  void FeedCentrals(HaloSnapshot_t &halo_snap);
  void BuildHDFDataType();
  void BuildMPIDataType();
public:
  SubhaloList_t Subhalos;
  MemberShipTable_t MemberTable;
  
  SubhaloSnapshot_t(): Snapshot_t(), Subhalos(), MemberTable(), ParallelizeHaloes(true)
  {
	BuildHDFDataType();
	BuildMPIDataType();
  }
  ~SubhaloSnapshot_t()
  {
	H5Tclose(H5T_SubhaloInDisk);
	H5Tclose(H5T_SubhaloInMem);
	MPI_Type_free(&MPI_HBT_SubhaloShell_t);
  }
  void GetSubFileName(string &filename, int iFile);
  void GetSrcFileName(string &filename, int iFile);
  void Load(mpi::communicator &world, int snapshot_index, bool load_src=false);
  void Save(mpi::communicator &world);
  void Clear()
  {
	//TODO
	cout<<"Clean() not implemented yet\n";
  }
  void UpdateParticles(mpi::communicator & world, const ParticleSnapshot_t & snapshot);
//   void ParticleIndexToId();
  void AverageCoordinates();
  void AssignHosts(mpi::communicator &world, HaloSnapshot_t &halo_snap, const ParticleSnapshot_t &part_snap);
  void PrepareCentrals(HaloSnapshot_t &halo_snap);
  void RefineParticles();
  void UpdateTracks(mpi::communicator &world, const HaloSnapshot_t &halo_snap);
  HBTInt size() const
  {
	return Subhalos.size();
  }
  HBTInt GetId(HBTInt index) const
  {
	return Subhalos[index].TrackId;
  }
  const HBTxyz & GetComovingPosition(HBTInt index) const
  {
	return Subhalos[index].ComovingPosition;
  }
  const HBTxyz & GetPhysicalVelocity(HBTInt index) const
  {
	return Subhalos[index].PhysicalVelocity;
  }
  HBTReal GetMass(HBTInt index) const
  {
	return Subhalos[index].Particles.size();
  }
};

inline HBTInt GetCoreSize(HBTInt nbound)
{
  int coresize=nbound*HBTConfig.SubCoreSizeFactor;
  if(coresize<HBTConfig.SubCoreSizeMin) coresize=HBTConfig.SubCoreSizeMin;
  if(coresize>nbound) coresize=nbound;
  return coresize;
}
#endif
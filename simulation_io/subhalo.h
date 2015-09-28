#ifndef SUBHALO_HEADER_INCLUDED
#define SUBHALO_HEADER_INCLUDED

#include <iostream>
#include <new>

#include "../datatypes.h"
#include "snapshot_number.h"
#include "halo.h"

class TrackParticle_t
{
public:
  HBTInt TrackId;
  HBTInt TrackParticleId;
  HBTInt SnapshotIndexOfLastIsolation; //the last snapshot when it was a central, only considering past snapshots.
  HBTInt SnapshotIndexOfLastMaxMass; //the snapshot when it has the maximum subhalo mass, only considering past snapshots.
  HBTInt SnapshotIndexOfBirth;//when the subhalo first becomes resolved
  HBTInt SnapshotIndexOfDeath;//when the subhalo first becomes un-resolved.
  HBTInt LastMaxMass;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  TrackParticle_t()
  {
	TrackId=-1;
	TrackParticleId=HBTInt();
	SnapshotIndexOfLastIsolation=SpecialConst::NullSnapshotId;
	SnapshotIndexOfLastMaxMass=SpecialConst::NullSnapshotId;
	SnapshotIndexOfBirth=SpecialConst::NullSnapshotId;
	SnapshotIndexOfDeath=SpecialConst::NullSnapshotId;
  }
  void SetTrackParticle(const HBTInt particle_id)
  {
	TrackParticleId=particle_id;
  }
  void SetTrackId(const HBTInt track_id)
  {
	TrackId=track_id;
  }
};

class SubHalo_t: public TrackParticle_t
{
  typedef ShallowList_t <HBTInt> ParticleShallowList_t;
public:
  ParticleShallowList_t Particles;
  HBTInt Nbound, Nsrc;
  HBTInt HostHaloId;
  HBTReal RmaxComoving;
  HBTReal VmaxPhysical;
  HBTReal RPoissonComoving;
  SubHalo_t(): TrackParticle_t(), Particles(), Nbound(0), Nsrc(0)
  {
  }
  void unbind()
  {//TODO
  }
};

class SubHaloSnapshot_t: public SnapshotNumber_t
{  
  typedef DeepList_t <HBTInt> ParticleList_t;
  typedef DeepList_t <SubHalo_t> SubHaloList_t;
public:
  ParticleList_t AllParticles;
  SubHaloList_t SubHalos;
  SubHaloSnapshot_t(): SnapshotNumber_t(), SubHalos(), AllParticles()
  {
  }
  void Load(Parameter_t &param, int snapshot_index)
  {//TODO
	cout<<"SubHaloSnapshot_t::Load() not implemented yet\n";
	if(SnapshotIndex<HBTConfig.MinSnapshotIndex)
	{// LoadNull();
	}
	else
	{
	}
  }
  void Save(Parameter_t &param)
  {
	//TODO
	cout<<"Save() not implemted yet\n";
  }
  void Clear()
  {
	//TODO
	cout<<"Clean() not implemented yet\n";
  }
  void ParticleIdToIndex(Snapshot_t & snapshot)
  {
	for(HBTInt i=0;i<AllParticles.Size();i++)
	  AllParticles[i]=snapshot.GetParticleIndex(AllParticles[i]);
  }
  void ParticleIndexToId(Snapshot_t & snapshot)
  {
	for(HBTInt i=0;i<AllParticles.Size();i++)
	  AllParticles[i]=snapshot.GetParticleId(AllParticles[i]);
  }
  void descend_into(HaloSnapshot_t &halo_snap, Snapshot_t &part_snap);
  void decide_centrals();
  void refine_particles();
};


#endif
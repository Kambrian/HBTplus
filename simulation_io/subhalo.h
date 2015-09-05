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
  HBTInt TrackParticleReference;
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
	TrackParticleReference=HBTInt();
	SnapshotIndexOfLastIsolation=SpecialConst::NullSnapshotId;
	SnapshotIndexOfLastMaxMass=SpecialConst::NullSnapshotId;
	SnapshotIndexOfBirth=SpecialConst::NullSnapshotId;
	SnapshotIndexOfDeath=SpecialConst::NullSnapshotId;
  }
  void SetTrackParticle(const HBTInt reference)
  {
	TrackParticleReference=reference;
  }
  void SetTrackId(const HBTInt track_id)
  {
	TrackId=track_id;
  }
};

class SubHalo_t: public Halo_t, public TrackParticle_t
{
public:
  HBTInt HostHaloId;
  SubHalo_t(): Halo_t(),TrackParticle_t()
  {
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
};

#endif
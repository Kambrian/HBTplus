#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <iostream>
#include <new>

#include "datatypes.h"
#include "simulation_io/simulation_io.h"

class Halo_t
{
public:
  List_t <ParticleReference_t> ParticleReferenceList;
  Halo_t(): ParticleReferenceList()
  {
  }
  void Load();
  void Clear();
};
class TrackParticle_t
{
public:
  HBTInt TrackId;
  ParticleReference_t TrackParticleReference;
  HBTInt SnapshotIndexOfLatestInfall; //the last snapshot when it was a central, only considering past snapshots.
  HBTInt SnapshotIndexOfLatestMaxMass; //the snapshot when it has the maximum subhalo mass, only considering past snapshots.
  HBTInt SnapshotIndexOfBirth;//when the subhalo first becomes resolved
  HBTInt SnapshotIndexOfDeath;//when the subhalo first becomes un-resolved.
  TrackParticle_t()
  {
	TrackId=-1;
	TrackParticleReference=ParticleReference_t();
	SnapshotIndexOfLatestInfall=SpecialConst::NullSnapshotId;
	SnapshotIndexOfLatestMaxMass=SpecialConst::NullSnapshotId;
	SnapshotIndexOfBirth=SpecialConst::NullSnapshotId;
	SnapshotIndexOfDeath=SpecialConst::NullSnapshotId;
  }
  void SetTrackParticle(const ParticleReference_t & particle)
  {
	TrackParticleReference=particle;
  }
  void SetTrackId(const HBTInt track_id)
  {
	TrackId=track_id;
  }
  Particle_t &RetrieveTrackParticleFromSnapshot(Snapshot_t & snapshot)
  {
	return snapshot.ParticleList[TrackParticleReference];
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

#endif
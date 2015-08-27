#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <iostream>
#include <new>

#include "../datatypes.h"
#include "snapshot_number.h"

class HaloSnapshot_t: public SnapshotNumber_t
{
  template <class PIDtype_t>
  void LoadGroupV3(Parameter_t &param, PIDtype_t dummy);
  void GetFileNameFormat(Parameter_t &param, string &format, int &FileCounts, bool &IsSubFile, bool &NeedByteSwap);
public:
  HBTInt NumberOfHaloes;
  HBTInt NumberOfParticles;
  HBTInt * ParticleList;
  HBTInt * HaloLength;
  HBTInt * HaloOffset;
  HaloSnapshot_t()
  {
	ParticleList=NULL;
  }
  void Load(Parameter_t &param, int snapshot_index);
  void Clear();
};

class TrackParticle_t
{
public:
  HBTInt TrackId;
  HBTInt TrackParticleReference;
  HBTInt SnapshotIndexOfLatestInfall; //the last snapshot when it was a central, only considering past snapshots.
  HBTInt SnapshotIndexOfLatestMaxMass; //the snapshot when it has the maximum subhalo mass, only considering past snapshots.
  HBTInt SnapshotIndexOfBirth;//when the subhalo first becomes resolved
  HBTInt SnapshotIndexOfDeath;//when the subhalo first becomes un-resolved.
  TrackParticle_t()
  {
	TrackId=-1;
	TrackParticleReference=HBTInt();
	SnapshotIndexOfLatestInfall=SpecialConst::NullSnapshotId;
	SnapshotIndexOfLatestMaxMass=SpecialConst::NullSnapshotId;
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
class Halo_t
{
public:
  List_t <HBTInt> Particles;
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
#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <iostream>
#include <new>

#include "../datatypes.h"
#include "snapshot_number.h"
#include "snapshot.h"

class Halo_t
{  
public:
  typedef vector <HBTInt> ParticleList_t;
  ParticleList_t Particles;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
};

class HaloSnapshot_t: public SnapshotNumber_t
{  
  typedef vector <Halo_t> HaloList_t;
  template <class PIDtype_t>
  void LoadGroupV3(PIDtype_t dummy);
  void GetFileNameFormat(string &format, int &FileCounts, bool &IsSubFile, bool &NeedByteSwap);
public:
  const Snapshot_t * SnapshotPointer;
  HaloList_t Halos;
  HBTInt TotNumberOfParticles;
  HBTInt NumPartOfLargestHalo;
  HaloSnapshot_t(): SnapshotNumber_t(), Halos(), SnapshotPointer(nullptr), TotNumberOfParticles(0), NumPartOfLargestHalo(0)
  {
  }
  void Load(int snapshot_index);
  void Clear();
  void ParticleIdToIndex(const Snapshot_t & snapshot);
  void ParticleIndexToId();
  void AverageCoordinates();
};

#endif
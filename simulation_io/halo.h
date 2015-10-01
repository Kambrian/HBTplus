#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <iostream>
#include <new>

#include "../datatypes.h"
#include "snapshot_number.h"
#include "snapshot.h"

class Halo_t
{  
  typedef vector <HBTInt> ParticleList_t;
public:
  ParticleList_t Particles;
  HBTxyz CenterOfMassComoving;
  HBTxyz AverageVelocityPhysical;
};

class HaloSnapshot_t: public SnapshotNumber_t
{  
  typedef vector <Halo_t> HaloList_t;
  template <class PIDtype_t>
  void LoadGroupV3(Parameter_t &param, PIDtype_t dummy);
  void GetFileNameFormat(Parameter_t &param, string &format, int &FileCounts, bool &IsSubFile, bool &NeedByteSwap);
public:
  HaloList_t Halos;
  HBTInt TotNumberOfParticles;
  HaloSnapshot_t(): SnapshotNumber_t(), Halos()
  {
  }
  void Load(Parameter_t &param, int snapshot_index);
  void Clear();
  void ParticleIdToIndex(Snapshot_t & snapshot);
  void ParticleIndexToId(Snapshot_t & snapshot);
  void AverageCoordinates(Snapshot_t & snapshot);
};

#endif
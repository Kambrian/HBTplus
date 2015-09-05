#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <iostream>
#include <new>

#include "../datatypes.h"
#include "snapshot_number.h"

class Halo_t
{  
  typedef ShallowList_t <HBTInt> ParticleShallowList_t;
public:
  ParticleShallowList_t Particles;
};

class HaloSnapshot_t: public SnapshotNumber_t
{  
  typedef DeepList_t <HBTInt> ParticleList_t;
  typedef DeepList_t <Halo_t> HaloList_t;
  template <class PIDtype_t>
  void LoadGroupV3(Parameter_t &param, PIDtype_t dummy);
  void GetFileNameFormat(Parameter_t &param, string &format, int &FileCounts, bool &IsSubFile, bool &NeedByteSwap);
public:
  ParticleList_t AllParticles;
  HaloList_t Halos;
  HaloSnapshot_t(): SnapshotNumber_t(), Halos(), AllParticles()
  {
  }
  void Load(Parameter_t &param, int snapshot_index);
  void Clear();
};

#endif
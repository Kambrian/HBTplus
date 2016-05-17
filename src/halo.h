#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <iostream>
#include <new>

#include "datatypes.h"
#include "snapshot_number.h"
#include "snapshot.h"

class Halo_t
{  
public:
  typedef vector <HBTInt> ParticleList_t;
  ParticleList_t Particles;
  HBTxyz ComovingAveragePosition;
  HBTxyz PhysicalAverageVelocity;
};

class HaloSnapshot_t: public Snapshot_t
{  
  typedef vector <Halo_t> HaloList_t;
public:
  HaloList_t Halos;
  HBTInt TotNumberOfParticles;
  HBTInt NumPartOfLargestHalo;
  const ParticleSnapshot_t * ParticleSnapshot;
  
  HaloSnapshot_t(): Snapshot_t(), Halos(), ParticleSnapshot(nullptr), TotNumberOfParticles(0), NumPartOfLargestHalo(0)
  {
  }
  void Load(int snapshot_index);
  void Clear();
  void ParticleIdToIndex(const ParticleSnapshot_t & snapshot);
  void ParticleIndexToId();
  void AverageCoordinates();
  HBTInt size() const
  { 
	return Halos.size();
  }
  const HBTxyz & GetComovingPosition(const HBTInt index) const
  {
	return Halos[index].ComovingAveragePosition;
  }
  const HBTxyz & GetPhysicalVelocity(const HBTInt index) const
  {
	return Halos[index].PhysicalAverageVelocity;
  }
  HBTReal GetMass(const HBTInt index) const
  {
	return Halos[index].Particles.size();
  }
};

#endif
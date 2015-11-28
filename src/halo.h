#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <iostream>
#include <new>
#include <algorithm>

#include "datatypes.h"
#include "snapshot_number.h"
#include "snapshot.h"

class Halo_t
{
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
	ar & Particles;
	ar & HaloId;
	ar & ComovingPosition;
	ar & PhysicalVelocity;
  }  
public:
  typedef vector <Particle_t> ParticleList_t;
  ParticleList_t Particles;
  HBTInt HaloId;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  void Relocate(mpi::communicator &world, int root, ParticleSnapshot_t &snap);
};

class HaloSnapshot_t: public Snapshot_t
{  
  typedef vector <Halo_t> HaloList_t;
  void ExchangeGroups(mpi::communicator &world);
public:
  const ParticleSnapshot_t * ParticleSnapshot;
  HaloList_t Halos;
  HBTInt TotNumberOfParticles;
  HBTInt NumPartOfLargestHalo;
  
  HaloSnapshot_t(): Snapshot_t(), Halos(), ParticleSnapshot(nullptr), TotNumberOfParticles(0), NumPartOfLargestHalo(0)
  {
  }
  void Load(mpi::communicator & world, int snapshot_index);
  void Clear();
  void ParticleIdToIndex(const ParticleSnapshot_t & snapshot);
  void ParticleIndexToId();
  void AverageCoordinates();
  HBTInt size() const
  { 
	return Halos.size();
  }
  const HBTxyz & GetComovingPosition(HBTInt index) const
  {
	return Halos[index].ComovingPosition;
  }
  const HBTxyz & GetPhysicalVelocity(HBTInt index) const
  {
	return Halos[index].PhysicalVelocity;
  }
  HBTReal GetMass(HBTInt index) const
  {
	return Halos[index].Particles.size();
  }
};

#endif
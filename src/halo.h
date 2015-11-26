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
  HBTInt HaloId;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
};

struct GroupFileSize_t
{
  HBTInt NumberOfGroups;
  HBTInt NumberOfParticles;
};

class HaloSnapshot_t: public Snapshot_t
{  
  typedef vector <Halo_t> HaloList_t;
  void ReadFile(int iFile, GroupFileSize_t &filesize, int read_level, HBTInt start_particle=0, HBTInt end_particle=-1);
  template <class PIDtype_t>
  void ReadFileV2V3(int iFile, GroupFileSize_t &filesize, int read_level, HBTInt start_particle, HBTInt end_particle);
  int CountFiles();
  void GetFileNameFormat(string &format, int &FileCounts, bool &IsSubFile, bool &NeedByteSwap);
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
#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <climits>
#include <iostream>
#include <new>
#include <algorithm>
#include <numeric>

#include "datatypes.h"
#include "snapshot_number.h"
#include "snapshot.h"
#include "mpi_wrapper.h"

class Halo_t
{
public:
  typedef vector <Particle_t> ParticleList_t;
  ParticleList_t Particles;
  HBTInt HaloId;
  HBTxyz ComovingAveragePosition;
  HBTxyz PhysicalAverageVelocity;
  HBTReal Mass;
  void AverageCoordinates();
  /* deprecated; use move assignment instead; 
   * shall not define destructor in order for default move to be implemented by the compiler.
  void MoveTo(Halo_t & dest, bool MoveParticle=true)
  {
	dest.HaloId=HaloId;
	copyHBTxyz(dest.ComovingPosition, ComovingPosition);
	copyHBTxyz(dest.PhysicalVelocity, PhysicalVelocity);
	if(MoveParticle)
	  dest.Particles.swap(Particles);
  }
  */
  HBTInt KickNullParticles();
};
extern void create_MPI_Halo_Id_type(MPI_Datatype &MPI_HBTHalo_Id_t);

class HaloSnapshot_t: public Snapshot_t
{  
  typedef vector <Halo_t> HaloList_t;
  MPI_Datatype MPI_HBT_HaloId_t;//MPI datatype ignoring the particle list
  void BuildMPIDataType();
public:
  HaloList_t Halos;
  HBTInt TotNumberOfParticles;
  HBTInt NumPartOfLargestHalo;
  MappedIndexTable_t<HBTInt, HBTInt> ParticleHash;
  
  HaloSnapshot_t(): Snapshot_t(), Halos(), TotNumberOfParticles(0), NumPartOfLargestHalo(0)
  {
	BuildMPIDataType();
  }
  ~HaloSnapshot_t()
  {
// 	Clear();
	My_Type_free(&MPI_HBT_HaloId_t);
  }
  void Load(MpiWorker_t & world, int snapshot_index);
  void Clear();
  void UpdateParticles(MpiWorker_t & world, const ParticleSnapshot_t & snapshot);
//   void ParticleIndexToId();
  void FillParticleHash();
  void ClearParticleHash();
  HBTInt size() const
  { 
	return Halos.size();
  }
  HBTInt GetId(HBTInt index) const
  {
	return Halos[index].HaloId;
  }
  const HBTxyz & GetComovingPosition(HBTInt index) const
  {
	return Halos[index].ComovingAveragePosition;
  }
  const HBTxyz & GetPhysicalVelocity(HBTInt index) const
  {
	return Halos[index].PhysicalAverageVelocity;
  }
  HBTReal GetMass(HBTInt index) const
  {
	return Halos[index].Mass;
  }
};

#endif

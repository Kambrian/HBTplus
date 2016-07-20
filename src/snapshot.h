#ifndef SNAPSHOT_H_INCLUDED
#define SNAPSHOT_H_INCLUDED

#include <iostream>
#include <sstream>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "datatypes.h"
#include "mymath.h"
#include "config_parser.h"
#include "snapshot_number.h"
#include "hash.h"
#include "mpi_wrapper.h"

struct Cosmology_t
{
  HBTReal OmegaM0;
  HBTReal OmegaLambda0;
  HBTReal ScaleFactor;
  
  //derived parameters:
  HBTReal Hz; //current Hubble param in internal units
  HBTReal OmegaZ;
  
  void Set(double scalefactor, double omega0, double omegaLambda0)
  {
	OmegaM0=omega0;
	OmegaLambda0=omegaLambda0;
	ScaleFactor=scalefactor;
	Hz=PhysicalConst::H0 * sqrt(OmegaM0 / (ScaleFactor * ScaleFactor * ScaleFactor) 
		+ (1 - OmegaM0 - OmegaLambda0) / (ScaleFactor * ScaleFactor)
		+ OmegaLambda0);//Hubble param for the current catalogue;
	
	HBTReal Hratio=Hz/PhysicalConst::H0;
	OmegaZ=OmegaM0/(ScaleFactor*ScaleFactor*ScaleFactor)/Hratio/Hratio;
  }
};

struct RadVelMass_t
{
  HBTReal r, v, m; 
};
  
struct Particle_t
{
  HBTInt Id;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  HBTReal Mass;
#ifndef DM_ONLY
#ifdef UNBIND_WITH_THERMAL_ENERGY
  HBTReal InternalEnergy;
#endif
  ParticleType_t Type;
#endif
  void create_MPI_type(MPI_Datatype &dtype);
  Particle_t()=default;
  Particle_t(HBTInt id): Id(id)
  {
  }
};
extern ostream& operator << (ostream& o, Particle_t &p);


class Snapshot_t: public SnapshotNumber_t
{
public:
  Cosmology_t Cosmology;
//   Snapshot_t()=default;
  virtual HBTInt size() const=0;
  virtual HBTInt GetId(const HBTInt index) const
  {
	return index;
  }
  virtual const HBTxyz & GetComovingPosition(const HBTInt index) const=0;
  virtual const HBTxyz & GetPhysicalVelocity(const HBTInt index) const=0;
  virtual HBTReal GetMass(const HBTInt index) const=0;
  virtual HBTReal GetInternalEnergy(HBTInt index) const
  {
	return 0.;
  }
  void SphericalOverdensitySize(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <HBTReal> &RSorted, HBTReal ParticleMass) const;
  void SphericalOverdensitySize(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <RadVelMass_t> &prof) const;
  void SphericalOverdensitySize2(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <HBTReal> &RSorted, HBTReal ParticleMass) const;
  void HaloVirialFactors(HBTReal &virialF_tophat, HBTReal &virialF_b200, HBTReal &virialF_c200) const;
  void RelativeVelocity(const HBTxyz& targetPos, const HBTxyz& targetVel, const HBTxyz& refPos, const HBTxyz& refVel, HBTxyz& relativeVel) const;
};

inline void Snapshot_t::RelativeVelocity(const HBTxyz& targetPos, const HBTxyz& targetVel, const HBTxyz& refPos, const HBTxyz& refVel, HBTxyz& relativeVel) const
{
   HBTxyz dx;
   HBTxyz &dv=relativeVel;
  for(int j=0;j<3;j++)
  {
	dx[j]=targetPos[j]-refPos[j];
	if(HBTConfig.PeriodicBoundaryOn)  dx[j]=NEAREST(dx[j]);
	dv[j]=targetVel[j]-refVel[j];
	dv[j]+=Cosmology.Hz*Cosmology.ScaleFactor*dx[j];
  }
}

class SnapshotView_t: public Snapshot_t
{
public:
  HBTInt * Ids;
  HBTInt N;
  Snapshot_t & Snapshot;
  SnapshotView_t(vector <HBTInt> & ids, Snapshot_t & fullsnapshot): Ids(ids.data()), N(ids.size()), Snapshot(fullsnapshot), Snapshot_t(fullsnapshot)
  {
  };
  SnapshotView_t(VectorView_t <HBTInt> &ids, Snapshot_t & fullsnapshot): Ids(ids.data()), N(ids.size()), Snapshot(fullsnapshot), Snapshot_t(fullsnapshot)
  {
  };
  SnapshotView_t(HBTInt *ids, HBTInt n, Snapshot_t & fullsnapshot): Ids(ids), N(n), Snapshot(fullsnapshot), Snapshot_t(fullsnapshot)
  {
  };
  void ReSize(HBTInt n)
  {
	N=n;
  }
  HBTInt size() const
  {
	return N;
  }
  HBTInt GetId(HBTInt i) const
  {
	return Snapshot.GetId(Ids[i]);
  }
  HBTReal GetMass(HBTInt i) const
  {
	return Snapshot.GetMass(Ids[i]);
  }
  const HBTxyz & GetPhysicalVelocity(HBTInt i) const
  {
	return Snapshot.GetPhysicalVelocity(Ids[i]);
  }
  const HBTxyz & GetComovingPosition(HBTInt i) const
  {
	return Snapshot.GetComovingPosition(Ids[i]);
  }
};

class ParticleSnapshot_t: public Snapshot_t
{
  typedef vector <HBTInt> IndexList_t;
   
  FlatIndexTable_t<HBTInt, HBTInt> FlatHash;
  MappedIndexTable_t<HBTInt, HBTInt> MappedHash;
  IndexTable_t<HBTInt, HBTInt> *ParticleHash;
  
  void ExchangeParticles(MpiWorker_t &world);
  void PartitionParticles(MpiWorker_t &world, vector <int> &offset);
  bool IsContiguousId(MpiWorker_t &world, HBTInt &GlobalIdMin);
  HBTInt IdMin, IdMax;
public:
  vector <Particle_t> Particles;
  HBTInt NumberOfParticlesOnAllNodes;
  vector <HBTInt> ProcessIdRanges;
  
  ParticleSnapshot_t(): Snapshot_t(), Particles(), ParticleHash(), MappedHash(), FlatHash(), NumberOfParticlesOnAllNodes(0)
  {
	if(HBTConfig.ParticleIdNeedHash)
	  ParticleHash=&MappedHash;
	else
	  ParticleHash=&FlatHash;
  }
  ParticleSnapshot_t(MpiWorker_t &world, int snapshot_index, bool fill_particle_hash=true): ParticleSnapshot_t()
  {
	Load(world, snapshot_index, fill_particle_hash);
  }
  ~ParticleSnapshot_t()
  {
	Clear();//not necessary
  }
  void FillParticleHash();
  void ClearParticleHash();
  
  HBTInt size() const;
  HBTInt GetId(HBTInt index) const;
  HBTInt GetIndex(HBTInt particle_id) const;
  HBTInt GetIndex(Particle_t & particle) const;
  template <class ParticleIdList_t>
  void GetIndices(ParticleIdList_t &particles) const;
  const HBTxyz & GetComovingPosition(HBTInt index) const;
  const HBTxyz & GetPhysicalVelocity(HBTInt index) const;
  HBTReal GetMass(HBTInt index) const;
  HBTReal GetInternalEnergy(HBTInt index) const;
  ParticleType_t GetParticleType(HBTInt index) const;
  
  void Load(MpiWorker_t &world, int snapshot_index, bool fill_particle_hash=true);
  void Clear();
  
  void AveragePosition(HBTxyz & CoM, const HBTInt Particles[], HBTInt NumPart) const; 
  void AverageVelocity(HBTxyz & CoV, const HBTInt Particles[], HBTInt NumPart) const;
  
  template <class Halo_T>
  void ExchangeHalos(MpiWorker_t &world, vector <Halo_T> & InHalos, vector <Halo_T> & OutHalos, MPI_Datatype MPI_Halo_Shell_Type) const;
};
inline HBTInt ParticleSnapshot_t::size() const
{
  return Particles.size();
}
inline HBTInt ParticleSnapshot_t::GetId(HBTInt index) const
{
  return Particles[index].Id;
}
inline HBTInt ParticleSnapshot_t::GetIndex(HBTInt particle_id) const
{
  return ParticleHash->GetIndex(particle_id);
}
inline HBTInt ParticleSnapshot_t::GetIndex(Particle_t & particle) const
{
  return ParticleHash->GetIndex(particle.Id);
}
inline const HBTxyz& ParticleSnapshot_t::GetComovingPosition(HBTInt index) const
{
  return Particles[index].ComovingPosition;
}
inline const HBTxyz& ParticleSnapshot_t::GetPhysicalVelocity(HBTInt index) const
{
  return Particles[index].PhysicalVelocity;
}
inline HBTReal ParticleSnapshot_t::GetMass(HBTInt index) const
{
  return Particles[index].Mass;
}
inline HBTReal ParticleSnapshot_t::GetInternalEnergy(HBTInt index) const
{
#if !defined(DM_ONLY) && defined(UNBIND_WITH_THERMAL_ENERGY) 
  return Particles[index].InternalEnergy;
#else
  return 0.;
#endif
}
inline ParticleType_t ParticleSnapshot_t::GetParticleType(HBTInt index) const
{
#ifdef DM_ONLY
  return TypeDM;
#else
  return Particles[index].Type;
#endif  
}

extern double AveragePosition(HBTxyz& CoM, const Particle_t Particles[], HBTInt NumPart);
extern double AverageVelocity(HBTxyz& CoV, const Particle_t Particles[], HBTInt NumPart);
#endif
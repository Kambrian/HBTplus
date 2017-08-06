#ifndef SNAPSHOT_H_INCLUDED
#define SNAPSHOT_H_INCLUDED

#include <iostream>
#include <sstream>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>

#include "datatypes.h"
#include "mymath.h"
#include "config_parser.h"
#include "snapshot_number.h"
#include "hash.h"

// #define TypeMax 6

struct RadVelMass_t
{
  HBTReal r, v, m; 
};
  
struct Particle_t
{
  HBTInt Id;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
#ifndef DM_ONLY
  HBTReal Mass;
#ifdef HAS_THERMAL_ENERGY
  HBTReal InternalEnergy;
#endif
  ParticleType_t Type;
#endif
  Particle_t(){};//do nothing. this leaves the content uninitialized, for fast memory allocation.
};

struct Cosmology_t
{
  HBTReal OmegaM0;
  HBTReal OmegaLambda0;
  HBTReal ScaleFactor;
  HBTReal ParticleMass;//DM particle mass if available.
  
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
  void HaloVirialFactors(HBTReal &virialF_tophat, HBTReal &virialF_b200, HBTReal &virialF_c200) const;
  void SphericalOverdensitySize(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <HBTReal> &RSorted) const;
  void SphericalOverdensitySize(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <RadVelMass_t> &prof) const;
  void SphericalOverdensitySize2(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <HBTReal> &RSorted) const;
};

class Snapshot_t: public SnapshotNumber_t
{
public:
  Cosmology_t Cosmology;
//   Snapshot_t()=default;
  virtual HBTInt size() const=0;
  virtual HBTInt GetMemberId(const HBTInt index) const
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
  void RelativeVelocity(const HBTxyz& targetPos, const HBTxyz& targetVel, const HBTxyz& refPos, const HBTxyz& refVel, HBTxyz& relativeVel) const;
};

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
  HBTInt ReSize(HBTInt n)
  {
	N=n;
  }
  HBTInt size() const
  {
	return N;
  }
  HBTInt GetMemberId(const HBTInt i) const
  {
	return Ids[i];
  }
  HBTReal GetMass(const HBTInt i) const
  {
	return Snapshot.GetMass(Ids[i]);
  }
  const HBTxyz & GetPhysicalVelocity(const HBTInt i) const
  {
	return Snapshot.GetPhysicalVelocity(Ids[i]);
  }
  const HBTxyz & GetComovingPosition(const HBTInt i) const
  {
	return Snapshot.GetComovingPosition(Ids[i]);
  }
};

class ParticleSnapshot_t: public Snapshot_t
{
  typedef HBTInt ParticleIndex_t ;
  typedef HBTInt ParticleId_t;
  typedef vector <ParticleIndex_t> IndexList_t;
  
  vector <Particle_t> Particles;
  FlatIndexTable_t<ParticleId_t, ParticleIndex_t> FlatHash;
  MappedIndexTable_t<ParticleId_t, ParticleIndex_t> MappedHash;
  IndexTable_t<ParticleId_t, ParticleIndex_t> *ParticleHash;
      
public:
  ParticleSnapshot_t(): Snapshot_t(), FlatHash(), MappedHash(), ParticleHash(), Particles()
  {	
	if(HBTConfig.ParticleIdNeedHash)
	  ParticleHash=&MappedHash;
	else
	  ParticleHash=&FlatHash;
  }
  ParticleSnapshot_t(int snapshot_index, bool fill_particle_hash=true): ParticleSnapshot_t()
  {
	Load(snapshot_index, fill_particle_hash);
  }
  ~ParticleSnapshot_t()
  {
	Clear();
  }
  HBTInt size() const;
  HBTInt GetMemberId(const HBTInt index) const;
  void FillParticleHash();
  void ClearParticleHash();
  ParticleIndex_t GetParticleIndex(const ParticleId_t particle_id) const;
  void Clear();
  ParticleIndex_t GetNumberOfParticles() const;
  ParticleId_t GetParticleId(const ParticleIndex_t index) const;
  const HBTxyz &GetComovingPosition(const ParticleIndex_t index) const;
  const HBTxyz &GetPhysicalVelocity(const ParticleIndex_t index) const;
  HBTReal GetParticleMass(const ParticleIndex_t index) const;
  HBTReal GetMass(const ParticleIndex_t index) const;
  HBTReal GetInternalEnergy(ParticleIndex_t index) const;
  ParticleType_t GetParticleType(ParticleIndex_t index) const;
  double AveragePosition(HBTxyz & CoM, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const; 
  double AverageVelocity(HBTxyz & CoV, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const;
  const Particle_t & GetParticle(HBTInt index) const;
  
  void Load(int snapshot_index, bool fill_particle_hash=true);
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

inline HBTInt ParticleSnapshot_t::size() const
{
  return Particles.size();
}
inline HBTInt ParticleSnapshot_t::GetMemberId(const HBTInt index) const
{
  return Particles[index].Id;
}
inline ParticleSnapshot_t::ParticleIndex_t ParticleSnapshot_t::GetParticleIndex(const ParticleId_t particle_id) const
{
//   return ParticleHash.at(particle_id);//this is safe
  return ParticleHash->GetIndex(particle_id);
}
inline ParticleSnapshot_t::ParticleIndex_t ParticleSnapshot_t::GetNumberOfParticles() const
{
  return size();
}
inline ParticleSnapshot_t::ParticleId_t ParticleSnapshot_t::GetParticleId(const ParticleIndex_t index) const
{
  return Particles[index].Id;
}
inline const HBTxyz& ParticleSnapshot_t::GetComovingPosition(const ParticleIndex_t index) const
{
  return Particles[index].ComovingPosition;
}
inline const HBTxyz& ParticleSnapshot_t::GetPhysicalVelocity(const ParticleIndex_t index) const
{
  return Particles[index].PhysicalVelocity;
}
inline HBTReal ParticleSnapshot_t::GetParticleMass(const ParticleIndex_t index) const
{
#ifdef DM_ONLY
	return Cosmology.ParticleMass;
#else
	return Particles[index].Mass;
#endif
}
inline HBTReal ParticleSnapshot_t::GetMass(const HBTInt index) const
{
  return GetParticleMass(index);
}
inline HBTReal ParticleSnapshot_t::GetInternalEnergy(HBTInt index) const
{
#if !defined(DM_ONLY) && defined(HAS_THERMAL_ENERGY) 
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
inline const Particle_t& ParticleSnapshot_t::GetParticle(HBTInt index) const
{
  return Particles[index];
}
#endif
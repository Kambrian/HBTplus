#ifndef SNAPSHOT_H_INCLUDED
#define SNAPSHOT_H_INCLUDED

#include <iostream>
#include <sstream>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include "../datatypes.h"
#include "../mymath.h"
#include "../config_parser.h"
#include "snapshot_number.h"
#include "hash.h"

#define NUMBER_OF_PARTICLE_TYPES 6
#define SNAPSHOT_HEADER_SIZE 256

struct SnapshotHeader_t
{
  int      npart[NUMBER_OF_PARTICLE_TYPES];
  double   mass[NUMBER_OF_PARTICLE_TYPES];
  double   ScaleFactor;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned npartTotal[NUMBER_OF_PARTICLE_TYPES];  //differ from standard. to be able to hold large integers
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[SNAPSHOT_HEADER_SIZE- NUMBER_OF_PARTICLE_TYPES*4- NUMBER_OF_PARTICLE_TYPES*8- 2*8- 2*4- NUMBER_OF_PARTICLE_TYPES*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};
class Snapshot_t: public SnapshotNumber_t
{
public:
  double Hz; //current Hubble param in internal units
  double ScaleFactor;
  Snapshot_t(): Hz(0.), ScaleFactor(0.), SnapshotNumber_t()
  {
  }
  void SetEpoch(double scalefactor, double Omega0, double OmegaLambda)
  {
	ScaleFactor=scalefactor;
	Hz=PhysicalConst::H0 * sqrt(Omega0 / (ScaleFactor * ScaleFactor * ScaleFactor) 
  + (1 - Omega0 - OmegaLambda) / (ScaleFactor * ScaleFactor)
  + OmegaLambda);//Hubble param for the current catalogue;
  }
  void SetEpoch(const Snapshot_t & snap)
  {
	ScaleFactor=snap.ScaleFactor;
	Hz=snap.Hz;
  }
  virtual HBTInt GetSize() const=0;
  virtual HBTInt GetMemberId(const HBTInt index) const
  {
	return index;
  }
  virtual const HBTxyz & GetComovingPosition(const HBTInt index) const=0;
  virtual const HBTxyz & GetPhysicalVelocity(const HBTInt index) const=0;
  virtual HBTReal GetMass(const HBTInt index) const=0;
};
class ParticleSnapshot_t: public Snapshot_t
{
  typedef HBTInt ParticleIndex_t ;
  typedef HBTInt ParticleId_t;
  typedef vector <ParticleIndex_t> IndexList_t;
    
  bool PeriodicBox;
  bool NeedByteSwap;
  int IntTypeSize;
  int RealTypeSize;
  vector <ParticleIndex_t> NumberOfDMParticleInFiles;
  vector <ParticleIndex_t> OffsetOfDMParticleInFiles;
  
  ParticleIndex_t NumberOfParticles;
  ParticleId_t * ParticleId; //better hide this from the user!!!!! 
  HBTxyz * ComovingPosition;
  HBTxyz * PhysicalVelocity;
  HBTReal * ParticleMass;
  FlatIndexTable_t<ParticleId_t, ParticleIndex_t> FlatHash;
  MappedIndexTable_t<ParticleId_t, ParticleIndex_t> MappedHash;
  IndexTable_t<ParticleId_t, ParticleIndex_t> *ParticleHash;
//   unordered_map <ParticleId_t, ParticleIndex_t> ParticleHash;//TODO: optimize this;also use intel concurrent_unordered_map
  
  void LoadId(Parameter_t & param);
  void LoadPosition(Parameter_t & param);
  void LoadVelocity(Parameter_t & param);
  void LoadMass(Parameter_t & param);
  void LoadHeader(Parameter_t & param, int ifile=1);
  bool ReadFileHeader(FILE *fp, SnapshotHeader_t &header);
  ParticleIndex_t ReadNumberOfDMParticles(Parameter_t & param, int ifile);
  size_t SkipBlock(FILE *fp);
  size_t ReadBlock(FILE *fp, void *block, const size_t n_read, const size_t n_skip_before=0, const size_t n_skip_after=0);
  void * LoadBlock(Parameter_t &param, int block_id, size_t element_size, int dimension=1, bool is_massblock=false);
public:
  SnapshotHeader_t Header;
  ParticleSnapshot_t(): Snapshot_t(), Header()
  {
	PeriodicBox=true;
	NeedByteSwap=false;
	IntTypeSize=0;
	RealTypeSize=0;
	NumberOfParticles=0;
	ParticleId=NULL;
	ComovingPosition=NULL;
	PhysicalVelocity=NULL;
	ParticleMass=NULL;
  }
  ~ParticleSnapshot_t()
  {
	Clear();
  }
  HBTInt GetSize() const;
  HBTInt GetMemberId(const HBTInt index) const;
  void FillParticleHash();
  void ClearParticleHash();
  ParticleIndex_t GetParticleIndex(const ParticleId_t particle_id) const;
  void GetFileName(Parameter_t &param, int ifile, string &filename);
  void Clear();
  ParticleIndex_t GetNumberOfParticles() const;
  ParticleId_t GetParticleId(const ParticleIndex_t index) const;
  const HBTxyz &GetComovingPosition(const ParticleIndex_t index) const;
  const HBTxyz &GetPhysicalVelocity(const ParticleIndex_t index) const;
  HBTReal GetParticleMass(const ParticleIndex_t index) const;
  HBTReal GetMass(const ParticleIndex_t index) const;
  void Load(int snapshot_index, bool load_id=true, bool load_position=true, bool load_velocity=true, bool load_mass=true, bool fill_particle_hash=true);
  void AveragePosition(HBTxyz & CoM, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const; 
  void AverageVelocity(HBTxyz & CoV, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const;
};
inline HBTInt ParticleSnapshot_t::GetSize() const
{
  return NumberOfParticles;
}
inline HBTInt ParticleSnapshot_t::GetMemberId(const HBTInt index) const
{
  return ParticleId[index];
}
inline ParticleSnapshot_t::ParticleIndex_t ParticleSnapshot_t::GetParticleIndex(const ParticleId_t particle_id) const
{
//   return ParticleHash.at(particle_id);//this is safe
  return ParticleHash->GetIndex(particle_id);
}
inline ParticleSnapshot_t::ParticleIndex_t ParticleSnapshot_t::GetNumberOfParticles() const
{
  return NumberOfParticles;
}
inline ParticleSnapshot_t::ParticleId_t ParticleSnapshot_t::GetParticleId(const ParticleIndex_t index) const
{
  return ParticleId[index];
}
inline const HBTxyz& ParticleSnapshot_t::GetComovingPosition(const ParticleIndex_t index) const
{
  return ComovingPosition[index];
}
inline const HBTxyz& ParticleSnapshot_t::GetPhysicalVelocity(const ParticleIndex_t index) const
{
  return PhysicalVelocity[index];
}
inline HBTReal ParticleSnapshot_t::GetParticleMass(const ParticleIndex_t index) const
{
  if(Header.mass[1])
	return Header.mass[1];
  else
	return ParticleMass[index];
}
inline HBTReal ParticleSnapshot_t::GetMass(const HBTInt index) const
{
  return GetParticleMass(index);
}
#endif
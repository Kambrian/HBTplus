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

#define NUMBER_OF_PARTICLE_TYPES 6
#define SNAPSHOT_HEADER_SIZE 256

struct SnapshotHeader_t
{
  int      npart[NUMBER_OF_PARTICLE_TYPES];
  double   mass[NUMBER_OF_PARTICLE_TYPES];
  double   time;
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
  double Hz;//current Hubble param in internal units, BT extension
  char     fill[SNAPSHOT_HEADER_SIZE- NUMBER_OF_PARTICLE_TYPES*4- NUMBER_OF_PARTICLE_TYPES*8- 2*8- 2*4- NUMBER_OF_PARTICLE_TYPES*4- 2*4 - 4*8 -8];  /* fills to 256 Bytes */
};
class Snapshot_t: public SnapshotNumber_t
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
  unordered_map <ParticleId_t, ParticleIndex_t> ParticleHash;//TODO: optimize this
  
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
  Snapshot_t(): SnapshotNumber_t(), Header()
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
  ~Snapshot_t()
  {
	Clear();
  }
  void FillParticleHash();
  void ClearParticleHash();
  ParticleIndex_t GetParticleIndex(ParticleId_t particle_id);
  void GetFileName(Parameter_t &param, int ifile, string &filename);
  void Clear();
  ParticleIndex_t GetNumberOfParticles() const;
  ParticleId_t GetParticleId(ParticleIndex_t index) const;
  HBTxyz &GetComovingPosition(ParticleIndex_t index) const;
  HBTxyz &GetPhysicalVelocity(ParticleIndex_t index) const;
  HBTReal GetParticleMass(ParticleIndex_t index) const;
  void Load(Parameter_t & param, int snapshot_index, bool load_id=true, bool load_position=true, bool load_velocity=true, bool load_mass=true, bool fill_particle_hash=true);
  void AveragePosition(HBTxyz & CoM, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const; 
  void AverageVelocity(HBTxyz & CoV, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const;
};

inline Snapshot_t::ParticleIndex_t Snapshot_t::GetParticleIndex(ParticleId_t particle_id)
{
  return ParticleHash[particle_id];
}
inline Snapshot_t::ParticleIndex_t Snapshot_t::GetNumberOfParticles() const
{
  return NumberOfParticles;
}
inline Snapshot_t::ParticleId_t Snapshot_t::GetParticleId(ParticleIndex_t index) const
{
  return ParticleId[index];
}
inline HBTxyz& Snapshot_t::GetComovingPosition(ParticleIndex_t index) const
{
  return ComovingPosition[index];
}
inline HBTxyz& Snapshot_t::GetPhysicalVelocity(ParticleIndex_t index) const
{
  return PhysicalVelocity[index];
}
inline HBTReal Snapshot_t::GetParticleMass(ParticleIndex_t index) const
{
  if(Header.mass[1])
	return Header.mass[1];
  else
	return ParticleMass[index];
}
#endif
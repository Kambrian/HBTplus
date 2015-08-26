#ifndef SIMULATION_IO_H_INCLUDED
#define SIMULATION_IO_H_INCLUDED

#include <iostream>
#include <sstream>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include "../datatypes.h"
#include "../config_parser.h"

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
class Snapshot_t
{
  SnapshotHeader_t Header;
  
  bool NeedByteSwap;
  int IntTypeSize;
  int RealTypeSize;
  vector <HBTInt> NumberOfDMParticleInFiles;
  vector <HBTInt> OffsetOfDMParticleInFiles;
  
  int SnapshotIndex;
  int SnapshotId; //original number
  HBTInt NumberOfParticles;
  HBTInt * ParticleId; //better hide this from the user!!!!! 
  HBTxyz * ComovingPosition;
  HBTxyz * PhysicalVelocity;
  HBTReal * ParticleMass;
  unordered_map <HBTInt, HBTInt> ParticleHash;//TODO: optimize this
  
  void SetSnapshotIndex(Parameter_t &param, int snapshot_index);
  void LoadId(Parameter_t & param);
  void LoadPosition(Parameter_t & param);
  void LoadVelocity(Parameter_t & param);
  void LoadMass(Parameter_t & param);
  void LoadHeader(Parameter_t & param, int ifile);
  bool ReadFileHeader(FILE *fp, SnapshotHeader_t &header);
  HBTInt ReadNumberOfDMParticles(Parameter_t & param, int ifile);
  size_t SkipBlock(FILE *fp);
  size_t ReadBlock(FILE *fp, void *buf, const size_t n_read, const size_t n_skip_before, const size_t n_skip_after);
  void * LoadBlock(Parameter_t &param, int block_id, size_t element_size, int dimension, bool is_massblock);
public:
  Snapshot_t()
  {
	NeedByteSwap=false;
	IntTypeSize=0;
	RealTypeSize=0;
	SnapshotIndex=SpecialConst::NullSnapshotId;
	NumberOfParticles=0;
	ParticleId=NULL;
	ComovingPosition=NULL;
	PhysicalVelocity=NULL;
	ParticleMass=NULL;
  }
  void FillParticleHash();
  void ClearParticleHash();
  HBTInt GetParticleIndex(HBTInt particle_id);
  void GetFileName(Parameter_t &param, int ifile, string &filename);
  void FormatSnapshotId(std::stringstream &ss);
  void Clear();
  int GetSnapshotIndex();
  int GetSnapshotId();
  HBTInt GetNumberOfParticles();
  HBTInt GetParticleId(HBTInt index);
  HBTxyz &GetComovingPosition(HBTInt index);
  HBTxyz &GetPhysicalVelocity(HBTInt index);
  HBTReal GetParticleMass(HBTInt index);
  void Load(int snapshot_index, Parameter_t & param, bool load_id, bool load_position, bool load_velocity, bool load_mass, bool fill_particle_hash);
};
inline HBTInt Snapshot_t::GetParticleIndex(HBTInt particle_id)
{
  return ParticleHash[particle_id];
}
inline int Snapshot_t::GetSnapshotIndex()
{
  return SnapshotIndex;
}
inline int Snapshot_t::GetSnapshotId()
{
  return SnapshotId;
}
inline void Snapshot_t::FormatSnapshotId(stringstream& ss)
{
  ss << std::setw(3) << std::setfill('0') << SnapshotId;
}
inline void Snapshot_t::SetSnapshotIndex(Parameter_t & param, int snapshot_index)
{
  assert(snapshot_index>=param.MinSnapshotIndex&&snapshot_index<=param.MaxSnapshotIndex);
//   assert(SpecialConst::NullSnapshotId!=snapshot_index);
  SnapshotIndex=snapshot_index; 
  if(param.SnapshotIdList.empty())
	SnapshotId=SnapshotIndex;
  else
	SnapshotId=param.SnapshotIdList[SnapshotIndex];
}
inline HBTInt Snapshot_t::GetNumberOfParticles()
{
  return NumberOfParticles;
}
inline HBTInt Snapshot_t::GetParticleId(HBTInt index)
{
  return ParticleId[index];
}
inline HBTxyz& Snapshot_t::GetComovingPosition(HBTInt index)
{
  return ComovingPosition[index];
}
inline HBTxyz& Snapshot_t::GetPhysicalVelocity(HBTInt index)
{
  return PhysicalVelocity[index];
}
inline HBTReal Snapshot_t::GetParticleMass(HBTInt index)
{
  if(Header.mass[1])
	return Header.mass[1];
  else
	return ParticleMass[index];
}
#endif
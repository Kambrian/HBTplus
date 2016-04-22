#ifndef SNAPSHOT_H_INCLUDED
#define SNAPSHOT_H_INCLUDED

#include <iostream>
#include <sstream>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include <list>
#include <forward_list>

#include "datatypes.h"
#include "mymath.h"
#include "config_parser.h"
#include "snapshot_number.h"
#include "hash.h"
#include "mpi_wrapper.h"

#define NUMBER_OF_PARTICLE_TYPES 6
#define SNAPSHOT_HEADER_SIZE 256
class Particle_t
{
public:
  HBTInt Id;
//   HBTInt ParticleIndex;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  HBTReal Mass;
  Particle_t()=default;
  Particle_t(HBTInt id): Id(id)
  {
  }
  void create_MPI_type(MPI_Datatype &MPI_HBTParticle_t);
};
extern ostream& operator << (ostream& o, Particle_t &p);

class SnapshotHeader_t
{
public:
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
  double   OmegaM0;
  double   OmegaLambda0;
  double   HubbleParam; 
  char     fill[SNAPSHOT_HEADER_SIZE- NUMBER_OF_PARTICLE_TYPES*4- NUMBER_OF_PARTICLE_TYPES*8- 2*8- 2*4- NUMBER_OF_PARTICLE_TYPES*4- 2*4 - 4*8];  /* fills to 256 Bytes */
  void create_MPI_type(MPI_Datatype &dtype);
};

class Snapshot_t: public SnapshotNumber_t
{
public:
  /*epoch header*/
  HBTReal OmegaM0;
  HBTReal OmegaLambda0;
  HBTReal Hz; //current Hubble param in internal units
  HBTReal ScaleFactor;
  /*end epoch header*/
  
  Snapshot_t(): Hz(0.), ScaleFactor(0.), OmegaM0(0.), OmegaLambda0(0.), SnapshotNumber_t()
  {
  }
  Snapshot_t(const Snapshot_t & sn)
  {
	SetEpoch(sn);
  }
  void SetEpoch(HBTReal scalefactor, HBTReal omegaM0, HBTReal omegaLambda0)
  {
	ScaleFactor=scalefactor;
	OmegaM0=omegaM0;
	OmegaLambda0=omegaLambda0;
	Hz=PhysicalConst::H0 * sqrt(OmegaM0 / (ScaleFactor * ScaleFactor * ScaleFactor) 
	+ (1 - OmegaM0 - OmegaLambda0) / (ScaleFactor * ScaleFactor)
	+ OmegaLambda0);//Hubble param for the current catalogue;
  }
  void SetEpoch(const Snapshot_t & snap)
  {
	ScaleFactor=snap.ScaleFactor;
	Hz=snap.Hz;
	OmegaM0=snap.OmegaM0;
	OmegaLambda0=snap.OmegaLambda0;
  }
  virtual HBTInt size() const=0;
  virtual HBTInt GetId(HBTInt index) const
  {
	return index;
  }
  virtual const HBTxyz & GetComovingPosition(HBTInt index) const=0;
  virtual const HBTxyz & GetPhysicalVelocity(HBTInt index) const=0;
  virtual HBTReal GetMass(HBTInt index) const=0;
  void SphericalOverdensitySize(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <HBTReal> &RSorted, HBTReal ParticleMass) const;
  void SphericalOverdensitySize2(float &Mvir, float &Rvir, HBTReal VirialFactor, const vector <HBTReal> &RSorted, HBTReal ParticleMass) const;
  void HaloVirialFactors(HBTReal &virialF_tophat, HBTReal &virialF_b200, HBTReal &virialF_c200) const;
};
class SnapshotView_t: public Snapshot_t
{
public:
  HBTInt * Ids;
  HBTInt N;
  Snapshot_t & Snapshot;
  SnapshotView_t(vector <HBTInt> & ids, Snapshot_t & fullsnapshot): Ids(ids.data()), N(ids.size()), Snapshot(fullsnapshot)
  {
	SetEpoch(fullsnapshot);
  };
  SnapshotView_t(VectorView_t <HBTInt> &ids, Snapshot_t & fullsnapshot): Ids(ids.data()), N(ids.size()), Snapshot(fullsnapshot)
  {
	SetEpoch(fullsnapshot);
  };
  SnapshotView_t(HBTInt *ids, HBTInt n, Snapshot_t & fullsnapshot): Ids(ids), N(n), Snapshot(fullsnapshot)
  {
	SetEpoch(fullsnapshot);
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
  /*header extension*/
  bool NeedByteSwap;
  int IntTypeSize;
  int RealTypeSize;
  vector <HBTInt> NumberOfDMParticleInFiles;
  vector <HBTInt> OffsetOfDMParticleInFiles;
    
  HBTInt NumberOfParticles;
  FlatIndexTable_t<HBTInt, HBTInt> FlatHash;
  MappedIndexTable_t<HBTInt, HBTInt> MappedHash;
  IndexTable_t<HBTInt, HBTInt> *ParticleHash;
  
  void ReadFile(int ifile);
  void LoadHeader(int ifile=0);
  bool ReadFileHeader(FILE *fp, SnapshotHeader_t &header);
  HBTInt ReadNumberOfDMParticles(int ifile);
  size_t SkipBlock(FILE *fp);
  void ExchangeParticles(MpiWorker_t &world);
public:
  SnapshotHeader_t Header;
  vector <Particle_t> Particles;
  
  ParticleSnapshot_t(): Snapshot_t(), Header(), Particles(), NumberOfParticles(0), ParticleHash(), MappedHash(), FlatHash()
  {
	NeedByteSwap=false;
	IntTypeSize=0;
	RealTypeSize=0;
	if(HBTConfig.ParticleIdNeedHash)
	  ParticleHash=&MappedHash;
	else
	  ParticleHash=&FlatHash;
  }
  ~ParticleSnapshot_t()
  {
	Clear();//not necessary
  }
  void FillParticleHash();
  void ClearParticleHash();
  void GetFileName(int ifile, string &filename);
  
  HBTInt size() const;
  HBTInt GetId(HBTInt index) const;
  HBTInt GetIndex(HBTInt particle_id) const;
  HBTInt GetIndex(Particle_t & particle) const;
  const HBTxyz & GetComovingPosition(HBTInt index) const;
  const HBTxyz & GetPhysicalVelocity(HBTInt index) const;
  HBTReal GetMass(HBTInt index) const;
  
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
{/*
  if(Header.mass[1])
	return Header.mass[1];
  else*/
	return Particles[index].Mass;
}

extern void AveragePosition(HBTxyz& CoM, const Particle_t Particles[], HBTInt NumPart);
extern void AverageVelocity(HBTxyz& CoV, const Particle_t Particles[], HBTInt NumPart);
#endif
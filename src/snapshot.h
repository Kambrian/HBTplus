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
#include "boost_mpi.h"

#define NUMBER_OF_PARTICLE_TYPES 6
#define SNAPSHOT_HEADER_SIZE 256

struct Particle_t
{
  HBTInt Id;
//   HBTInt ParticleIndex;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  HBTReal Mass;
};

class SnapshotHeader_t
{
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);//for boost 
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
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[SNAPSHOT_HEADER_SIZE- NUMBER_OF_PARTICLE_TYPES*4- NUMBER_OF_PARTICLE_TYPES*8- 2*8- 2*4- NUMBER_OF_PARTICLE_TYPES*4- 2*4 - 4*8];  /* fills to 256 Bytes */
  
};
template<class Archive>
void SnapshotHeader_t::serialize(Archive& ar, const unsigned int version)
{
    ar &  npart;
    ar &  mass;
    ar &  ScaleFactor;
    ar &  redshift;
    ar &  flag_sfr;
    ar &  flag_feedback;
    ar &  npartTotal;
    ar &  flag_cooling;
    ar &  num_files;
    ar &  BoxSize;
    ar &  Omega0;
    ar &  OmegaLambda;
    ar &  HubbleParam;
    ar &  fill; 
}
BOOST_IS_MPI_DATATYPE(SnapshotHeader_t)

class Snapshot_t: public SnapshotNumber_t
{
public:
  /*epoch header*/
  double Hz; //current Hubble param in internal units
  double ScaleFactor;
  /*end epoch header*/
  
  Snapshot_t(): Hz(0.), ScaleFactor(0.), SnapshotNumber_t()
  {
  }
  Snapshot_t(const Snapshot_t & sn)
  {
	SetEpoch(sn);
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
  virtual HBTInt size() const=0;
  virtual HBTInt GetId(HBTInt index) const
  {
	return index;
  }
  virtual const HBTxyz & GetComovingPosition(HBTInt index) const=0;
  virtual const HBTxyz & GetPhysicalVelocity(HBTInt index) const=0;
  virtual HBTReal GetMass(HBTInt index) const=0;
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
  vector <Particle_t> Particle;
  FlatIndexTable_t<HBTInt, HBTInt> FlatHash;
  MappedIndexTable_t<HBTInt, HBTInt> MappedHash;
  IndexTable_t<HBTInt, HBTInt> *ParticleHash;
  
  void ReadFile(int ifile);
  void LoadHeader(int ifile=0);
  bool ReadFileHeader(FILE *fp, SnapshotHeader_t &header);
  HBTInt ReadNumberOfDMParticles(int ifile);
  size_t SkipBlock(FILE *fp);
  void ExchangeParticles(mpi::communicator &world);
public:
  SnapshotHeader_t Header;
  ParticleSnapshot_t(): Snapshot_t(), Header(), Particle(), NumberOfParticles(0), ParticleHash()
  {
	NeedByteSwap=false;
	IntTypeSize=0;
	RealTypeSize=0;
  }
  ~ParticleSnapshot_t()
  {
	Clear();
  }
  void FillParticleHash();
  void ClearParticleHash();
  void GetFileName(int ifile, string &filename);
  
  HBTInt size() const;
  HBTInt GetId(HBTInt index) const;
  HBTInt GetIndex(HBTInt particle_id) const;
  const HBTxyz & GetComovingPosition(HBTInt index) const;
  const HBTxyz & GetPhysicalVelocity(HBTInt index) const;
  HBTReal GetMass(HBTInt index) const;
  
  void Load(mpi::communicator &world, int snapshot_index, bool fill_particle_hash=true);
  void Clear();
  
  void AveragePosition(HBTxyz & CoM, const HBTInt Particles[], HBTInt NumPart) const; 
  void AverageVelocity(HBTxyz & CoV, const HBTInt Particles[], HBTInt NumPart) const;
};
inline HBTInt ParticleSnapshot_t::size() const
{
  return Particle.size();
}
inline HBTInt ParticleSnapshot_t::GetId(HBTInt index) const
{
  return Particle[index].Id;
}
inline HBTInt ParticleSnapshot_t::GetIndex(HBTInt particle_id) const
{
  return ParticleHash->GetIndex(particle_id);
}
inline const HBTxyz& ParticleSnapshot_t::GetComovingPosition(HBTInt index) const
{
  return Particle[index].ComovingPosition;
}
inline const HBTxyz& ParticleSnapshot_t::GetPhysicalVelocity(HBTInt index) const
{
  return Particle[index].PhysicalVelocity;
}
inline HBTReal ParticleSnapshot_t::GetMass(HBTInt index) const
{
  if(Header.mass[1])
	return Header.mass[1];
  else
	return Particle[index].Mass;
}
#endif
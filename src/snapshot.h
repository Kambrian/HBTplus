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
  virtual HBTInt GetMemberId(const HBTInt index) const
  {
	return index;
  }
  virtual const HBTxyz & GetComovingPosition(const HBTInt index) const=0;
  virtual const HBTxyz & GetPhysicalVelocity(const HBTInt index) const=0;
  virtual HBTReal GetMass(const HBTInt index) const=0;
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
  /*header extension*/
  bool NeedByteSwap;
  int IntTypeSize;
  int RealTypeSize;
  vector <ParticleIndex_t> NumberOfDMParticleInFiles;
  vector <ParticleIndex_t> OffsetOfDMParticleInFiles;
  
  struct LoadFlag_t
  {
	bool Id;
	bool Pos;
	bool Vel;
	bool Mass;
	LoadFlag_t(): Id(true), Pos(true), Vel(true), Mass(true)
	{
	}
  } LoadFlag;
  
  ParticleIndex_t NumberOfParticles;
  vector <ParticleId_t> ParticleId; //better hide this from the user!!!!! 
  vector <HBTxyz> ComovingPosition;
  vector <HBTxyz> PhysicalVelocity;
  vector <HBTReal> ParticleMass;
  FlatIndexTable_t<ParticleId_t, ParticleIndex_t> FlatHash;
  MappedIndexTable_t<ParticleId_t, ParticleIndex_t> MappedHash;
  IndexTable_t<ParticleId_t, ParticleIndex_t> *ParticleHash;
  //   unordered_map <ParticleId_t, ParticleIndex_t> ParticleHash;//TODO: optimize this;also use intel concurrent_unordered_map
  
  void ReadFile(int ifile);
  template <class T, class U>
  void ReadScalarBlock(FILE *fp, size_t n_read, size_t n_skip, vector <U> &x);
  template <class T>
  void ReadXyzBlock(FILE *fp, size_t n_read, size_t n_skip, vector <HBTxyz> &x);
  void LoadPosition();
  void LoadVelocity();
  void LoadMass();
  void LoadHeader(int ifile=0);
  bool ReadFileHeader(FILE *fp, SnapshotHeader_t &header);
  ParticleIndex_t ReadNumberOfDMParticles(int ifile);
  size_t SkipBlock(FILE *fp);
  void * LoadBlock(int block_id, size_t element_size, bool is_massblock=false);
public:
  SnapshotHeader_t Header;
  ParticleSnapshot_t(): Snapshot_t(), Header(), LoadFlag(), ParticleId(), ComovingPosition(), PhysicalVelocity(), ParticleMass(), NumberOfParticles(0)
  {
	NeedByteSwap=false;
	IntTypeSize=0;
	RealTypeSize=0;
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
  void GetFileName(int ifile, string &filename);
  void Clear();
  ParticleIndex_t GetNumberOfParticles() const;
  ParticleId_t GetParticleId(const ParticleIndex_t index) const;
  const HBTxyz &GetComovingPosition(const ParticleIndex_t index) const;
  const HBTxyz &GetPhysicalVelocity(const ParticleIndex_t index) const;
  HBTReal GetParticleMass(const ParticleIndex_t index) const;
  HBTReal GetMass(const ParticleIndex_t index) const;
  void Load(int snapshot_index, bool fill_particle_hash=true);
  void SetLoadFlags(bool load_id, bool load_pos, bool load_vel, bool load_mass);
  void AveragePosition(HBTxyz & CoM, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const; 
  void AverageVelocity(HBTxyz & CoV, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const;
};
inline HBTInt ParticleSnapshot_t::size() const
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
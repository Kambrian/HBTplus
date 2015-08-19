#ifndef SIMULATION_IO_H_INCLUDED
#define SIMULATION_IO_H_INCLUDED

#include <iostream>
#include <sstream>
#include "../datatypes.h"

#define NumberOfParticleTypes 6

struct SnapshotHeader_t
{
  int      npart[NumberOfParticleTypes];
  double   mass[NumberOfParticleTypes];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned npartTotal[NumberOfParticleTypes];  //differ from standard. to be able to hold large integers
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  double Hz;//current Hubble param in internal units, BT extension
  int SnapshotIndex; //current snapshot
  char     fill[256- NumberOfParticleTypes*4- NumberOfParticleTypes*8- 2*8- 2*4- NumberOfParticleTypes*4- 2*4 - 4*8 -8 - 4];  /* fills to 256 Bytes */
};

class Snapshot_t
{
  SnapshotHeader_t Header;
  HBTInt NumberOfParticles;
  HBTInt * ParticleId; //better hide this from the user!!!!! 
  HBTxyz * ComovingPosition;
  HBTxyz * PhysicalVelocity;
  bool ByteOrder;
  void CheckSnapshotIndexIsValid();
  void LoadId();
  void LoadPosition();
  void LoadVelocity();
  void LoadHeader(int iFile);
public:
  Snapshot_t()
  {
	Header.SnapshotIndex=SpecialConst::NullSnapshotId;
	NumberOfParticles=0;
	ParticleId=nullptr;
	ComovingPosition=nullptr;
	PhysicalVelocity=nullptr;
  }
  void FormatSnapshotNumber(std::stringstream &ss);
  void Clear();
  void SetSnapshotIndex(int snapshot_index);
  int GetSnapshotIndex();
  HBTInt GetParticleId(HBTInt index);
  HBTxyz &GetComovingPosition(HBTInt index);
  HBTxyz &GetPhysicalVelocity(HBTInt index);
  HBTInt GetNumberOfParticles();
  void Load(bool load_id=true, bool load_position=true, bool load_velocity=true);
};

#endif
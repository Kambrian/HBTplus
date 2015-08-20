#ifndef SIMULATION_IO_H_INCLUDED
#define SIMULATION_IO_H_INCLUDED

#include <iostream>
#include <sstream>
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
  int SnapshotIndex;
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
  bool ReadFileHeader(FILE *fp, SnapshotHeader_t &header);
public:
  Snapshot_t()
  {
	SnapshotIndex=SpecialConst::NullSnapshotId;
	NumberOfParticles=0;
	ParticleId=NULL;
	ComovingPosition=NULL;
	PhysicalVelocity=NULL;
  }
  void GetSnapshotFileName(Parameter_t &param, int ifile, string &filename);
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
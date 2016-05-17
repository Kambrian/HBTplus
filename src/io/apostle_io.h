#ifndef APOSTLE_IO_INCLUDED
#define APOSTLE_IO_INCLUDED
#include "../hdf_wrapper.h"
#include "../halo.h"

struct ApostleHeader_t
{
  int      NumberOfFiles;
  double   BoxSize;
  double   ScaleFactor;
  double   OmegaM0;
  double   OmegaLambda0;
  double   mass[TypeMax];
  int      npart[TypeMax];  
  HBTInt npartTotal[TypeMax]; 
};

struct ParticleHost_t
{
  HBTInt ParticleId;
  HBTInt HostId;
};

class ApostleReader_t
{
  string SnapshotName;
    
  vector <HBTInt> np_file;
  vector <HBTInt> offset_file;
  ApostleHeader_t Header;
  void ReadHeader(int ifile, ApostleHeader_t &header);
  HBTInt CompileFileOffsets(int nfiles);
  void ReadSnapshot(int ifile, Particle_t * ParticlesInFile);
  void ReadGroupId(int ifile, ParticleHost_t * ParticlesInFile, bool FlagReadParticleId);
  void GetFileName(int ifile, string &filename);
  void SetSnapshot(int snapshotId);
  
  hid_t CountTable_t, MassTable_t, H5T_HBTxyz;
public:
  ApostleReader_t();
  void LoadSnapshot(int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology);
  HBTInt LoadGroups(int snapshotId, vector <Halo_t> &Halos);
  ~ApostleReader_t();
};

extern bool IsApostleGroup(const string &GroupFileFormat);
#endif
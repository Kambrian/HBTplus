/* IO for Apostle (EAGLE local group) data.
 * 
 * To specify a list of snapshot, list the snapshot directories (one per line) in snapshotlist.txt and place it under your subhalo output directory. 
 * 
 * To use this IO, in the config file, set SnapshotFormat to apostle,  and set GroupFileFormat to apostle or apostle_particle_index.
 * 
 * The groups loaded are already filled with particle properties, and the halos are distributed to processors according to the CoM of each halo.
 */

#ifndef APOSTLE_IO_INCLUDED
#define APOSTLE_IO_INCLUDED
#include "../hdf_wrapper.h"
#include "../halo.h"
#include "../mpi_wrapper.h"

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

void create_ApostleHeader_MPI_type(MPI_Datatype &dtype);

struct ParticleHost_t: public Particle_t
{
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
  void ReadGroupParticles(int ifile, ParticleHost_t * ParticlesInFile, bool FlagReadParticleId);
  void GetFileName(int ifile, string &filename);
  void SetSnapshot(int snapshotId);
  void GetParticleCountInFile(hid_t file, int np[]);
  void ExchangeAndMerge(MpiWorker_t &world, vector< Halo_t >& Halos);
  
  MPI_Datatype MPI_ApostleHeader_t;
  
public:
  ApostleReader_t()
  {
    create_ApostleHeader_MPI_type(MPI_ApostleHeader_t);
  }
  ~ApostleHeader_t()
  {
    MPI_Type_free(&MPI_ApostleHeader_t);
  }
  void LoadSnapshot(MpiWorker_t &world, int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology);
  HBTInt LoadGroups(MpiWorker_t &world, int snapshotId, vector <Halo_t> &Halos);
};

extern bool IsApostleGroup(const string &GroupFileFormat);
#endif
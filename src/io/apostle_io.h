#ifndef APOSTLE_IO_INCLUDED
#define APOSTLE_IO_INCLUDED
#include "../hdf_wrapper.h"

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

class ApostleReader_t
{
  int SnapshotId;
  string SnapshotName;
  vector <Particle_t> &Particles;
  Cosmology_t &Cosmology;
  
  vector <HBTInt> np_file;
  vector <HBTInt> offset_file;
  ApostleHeader_t Header;
  void ReadHeader(int ifile, ApostleHeader_t &header);
  HBTInt CompileFileOffsets(int nfiles);
  void LoadSnapshot();
  void ReadFile(int ifile);
  void GetFileName(int ifile, string &filename);
  
  hid_t CountTable_t, MassTable_t, H5T_HBTxyz;
public:
  ApostleReader_t(int snapshotId, vector <Particle_t> &particles, Cosmology_t &cosmology);
  ~ApostleReader_t();
};

#endif
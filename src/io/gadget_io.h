#ifndef GADGET_IO_INCLUDED
#define GADGET_IO_INCLUDED
#include "../mpi_wrapper.h"

#define SNAPSHOT_HEADER_SIZE 256

struct GadgetHeader_t
{
  int      npart[TypeMax];
  double   mass[TypeMax];
  double   ScaleFactor;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned npartTotal[TypeMax];  //differ from standard. to be able to hold large integers
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   OmegaM0;
  double   OmegaLambda0;
  double   HubbleParam; 
  char     fill[SNAPSHOT_HEADER_SIZE- TypeMax*4- TypeMax*8- 2*8- 2*4- TypeMax*4- 2*4 - 4*8];  /* fills to 256 Bytes */
  void create_MPI_type(MPI_Datatype& dtype);
};
  
class GadgetReader_t
{
  vector <Particle_t> &Particles;
  Cosmology_t &Cosmology;

  int SnapshotId;
  GadgetHeader_t Header;
  bool NeedByteSwap;
  int IntTypeSize;
  int RealTypeSize;
  vector <HBTInt> NumberOfParticleInFiles;
  vector <HBTInt> OffsetOfParticleInFiles;
  void ReadGadgetFile(int ifile);
  void LoadGadgetHeader(int ifile=0);
  bool ReadGadgetFileHeader(FILE *fp, GadgetHeader_t &header);
  HBTInt ReadGadgetNumberOfParticles(int ifile);
  void GetGadgetFileName(int ifile, string &filename);
  void Load(MpiWorker_t &world);
    
public:
  GadgetReader_t(MpiWorker_t &world, int snapshot_id, vector <Particle_t> &particles, Cosmology_t & cosmology);
};

#endif
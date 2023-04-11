/* IO for Gadget4 hdf data.
 *
 * Currently only designed/tested for Jiutian data which is DM-only;
 * support for multiple particle types to be added.
 *
 * To specify a list of snapshot, list the snapshot directories (one per line) in snapshotlist.txt and place it under your subhalo output directory.
 *
 * To use this IO, in the config file, set SnapshotFormat to gadget4_hdf,  and set GroupFileFormat to gadget4_hdf (or gadget4_hdf2, which is equivalent but uses a different algorithm for reading groups internally, and may be slower).
 *
 * The groups loaded are already filled with particle properties, and the halos are distributed to processors according to the CoM of each halo.
 */

#ifndef GADGET4_IO_INCLUDED
#define GADGET4_IO_INCLUDED
#include "../hdf_wrapper.h"
#include "../halo.h"
#include "../mpi_wrapper.h"

namespace Gadget4Reader
{
struct Gadget4Header_t
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

void create_Gadget4Header_MPI_type(MPI_Datatype &dtype);

struct ParticleHost_t: public Particle_t
{
  HBTInt HostId;
};

class Gadget4Reader_t
{
  string SnapshotName;

  vector <HBTInt> np_file;
  vector <HBTInt> offset_file;
  vector <HBTInt> nhalo_per_groupfile;
  vector <HBTInt> offsethalo_per_groupfile;

  const int root_node=0;
  //snap tab:
  vector <HBTInt> ProcLen;
  void CollectProcSizes(MpiWorker_t &world, const ParticleSnapshot_t &PartSnap);
  //group tab:
  vector <HBTInt> HaloSizesLocal;//halosizes read into this proc
  vector <HBTInt> HaloSizesAll;//only significant on root proc
  void LoadGroupTab(MpiWorker_t &world);

  Gadget4Header_t Header; //this is snapshot header, not necessarily filled when loading groups
  void ReadHeader(int ifile, Gadget4Header_t &header);
  int ReadGroupFileCounts(int ifile);
  HBTInt CompileFileOffsets(int nfiles);
  HBTInt CompileGroupFileOffsets(int nfiles);
  void ReadSnapshot(int ifile, Particle_t * ParticlesInFile);
  void ReadGroupLen(int ifile, HBTInt *buf);
  void GetFileName(int ifile, string &filename);
  void GetGroupFileName(int ifile, string &filename);
  void SetSnapshot(int snapshotId);
  void GetParticleCountInFile(hid_t file, int np[]);

  void LoadLeadingGroups(MpiWorker_t &world, const vector<Particle_t> &Particles, vector <Halo_t> &Halos);
  void LoadLocalGroups(MpiWorker_t &world, const vector<Particle_t> &Particles, vector <Halo_t> &Halos);

  MPI_Datatype MPI_Gadget4Header_t;

public:
  Gadget4Reader_t()
  {
    create_Gadget4Header_MPI_type(MPI_Gadget4Header_t);
  }
  ~Gadget4Reader_t()
  {
    My_Type_free(&MPI_Gadget4Header_t);
  }

  void LoadSnapshot(MpiWorker_t &world, int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology);
  void LoadGroups(MpiWorker_t &world, const ParticleSnapshot_t &partsnap, vector <Halo_t> &Halos);
};

extern bool IsGadget4Group(const string &GroupFileFormat);
}
#endif

/* IO for Apostle (EAGLE local group) data.
 *
 * To specify a list of snapshot, list the snapshot directories (one per line) in snapshotlist.txt and place it under your subhalo output directory.
 * 
 * This IO supports the following formats:
 *   SnapshotFormat: apostle, illustris
 *   GroupFileFormat: apostle, apostle_particle_index,  //for eagle or apostle data
 *                    illustris, illustris_particle_index //for illustris or illustrisTNG data.
 *                    apostle_helucid, apostle_helucid_particle_index//for helucid data. the halo_id of helucid starts from 1
 * 
 *          the *_particle_index format means the halos will be filled with particle indices in the snapshot, instead of particle ids, upon loading.
 * 
 * specifying `apostle_particle_index` or `illustris_particle_index` format will cause ApostleReader to fill the groups with particle indices that directly maps to the snapshot, saving any effort for translating from id to index later. If it is disirable to fill the groups with actual particle ids, please use `apostle` or `illustris` for the group format.
 *
 * The particles are loaded into snapshot memory type by type, loading first type file by file and then next type..
 */

#ifndef APOSTLE_IO_INCLUDED
#define APOSTLE_IO_INCLUDED
#include "../hdf_wrapper.h"
#include "../halo.h"

namespace Apostle{

typedef array <HBTInt, TypeMax> TypeCounts_t;

typedef enum 
{
    apostle,
    illustris,
    helucid
} SnapshotType_t;

struct ApostleSnap_t
{
  SnapshotType_t SnapType;
//   bool SnapIsIllustris;
//   int SnapshotId;
  string SnapDirBaseName;
  string SnapshotName;
  
  void SetSnapshot(int snapshotId);
};

class ApostleHeader_t: public ApostleSnap_t
{
  void GetParticleCountInFile(hid_t file, HBTInt np[]);
  void ReadFileHeader(int ifile);
  void CompileTypeOffsets();
public:    
  int      NumberOfFiles;
  double   BoxSize;
  double   ScaleFactor;
  double   OmegaM0;
  double   OmegaLambda0;
  double   Mass[TypeMax];
  TypeCounts_t  NumPart;
  TypeCounts_t NumPartTotal;
  HBTInt NumPartAll;
  
  HBTInt NullGroupId;
  HBTInt MinGroupId;
  
  TypeCounts_t TypeOffsetInMem; //offsets (in particles) of each type in the global particle array
  vector<TypeCounts_t> NumPartEachFile; //TypeCount in each file 
  vector<TypeCounts_t> FileOffsetInType; //cumsum(NumPartTypeFile, axis=file), sum over previous files for each type
  
  void GetFileName(int ifile, string &filename);
  void Fill(int snapshotId);
//   void LoadExtraHeaderParams();
};

class IllustrisGroupHeader_t: public ApostleSnap_t
{
    void ReadFileHeader(int ifile);
    void CompileFileOffsets();
public:
    int NumFiles;
    int NumGroupsThisFile;
    int NumGroupsTotal;
    
    vector<HBTInt> GroupFileLen, GroupFileOffset; //number of groups and offsets in groups for each file
    
    void GetFileName(int ifile, string &filename);
    void Fill(int snapshotId);
};

struct ParticleHost_t
{
  HBTInt ParticleId;
  HBTInt HostId;
};

class ApostleReader_t
{
//   const int NullGroupId=1<<30; //1073741824

  ApostleHeader_t SnapHeader;
  IllustrisGroupHeader_t GroupHeader;
  
  void ReadSnapshot(int ifile, Particle_t * Particles);
  void ReadGroupId(int ifile, ParticleHost_t *Particles, bool FlagReadParticleId);
  void ReadParticleIDs(int ifile, HBTInt *ParticleIDs);
public:
  void LoadSnapshot(int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology);
  HBTInt LoadApostleGroups(int snapshotId, vector <Halo_t> &Halos);
  HBTInt LoadIllustrisGroups(int snapshotId, vector <Halo_t> &Halos);
};

extern bool IsHelucidGroup(const string &GroupFileFormat);
extern bool IsHelucidSnap(const string &SnapshotFormat);
extern bool IsApostleGroup(const string &GroupFileFormat);
extern bool IsIllustrisGroup(const string &GroupFileFormat);
extern bool IsApostleSnap(const string &SnapshotFormat);
extern bool IsIllustrisSnap(const string &SnapshotFormat);
}
#endif

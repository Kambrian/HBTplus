#ifndef JING_IO_INCLUDED
#define JING_IO_INCLUDED

#include "../../halo.h"

struct JingHeader_t
{
   HBTInt Np;   //number of particles in the simulation
   HBTInt ips;  //time step number, or snapshot num
   HBTReal Redshift; // current redshift
   HBTReal Omegat;  // current omega_m, mass density
   HBTReal Lambdat; // current Omega_lambda, dark energy density
   HBTReal BoxSize;  // boxsize in unit of Mpc/h
   HBTReal xscale;
   HBTReal vscale;
   HBTReal Hz; //Hubble Param at ztp
   HBTReal vunit; //velocity unit to get physical peculiar velocity in km/s
   //extra parameters not written in snapshotfile
   HBTReal OmegaM0;  //omega_m0
   HBTReal OmegaLambda0; //OmegaLambda0;   
   HBTReal RedshiftIni;
   HBTInt SnapDivScale;
   bool FlagHasScale;
   bool ParticleDataXMajor;
   /*==extension needed for calculating binding energy:==*/
   HBTReal ScaleFactor;//current reduced scale factor 
   HBTReal mass[2];//Gas (mass[0]) and DM (mass[1]) particle masses, in units of 10^10Msun/h
};
class JingReader_t
{
  int SnapshotId;
  JingHeader_t Header;
  bool NeedByteSwap;
  int NumFilesPos, NumFilesVel, NumFilesId;
  
  string GetFileName(const char * filetype, int iFile=0);
  int CountFiles(const char *filetype); 
  void ProbeFiles();
  void ReadId(vector <Particle_t> &Particles);
  void ReadVelocity(vector <Particle_t> &Particles);
  void ReadPosition(vector <Particle_t> &Particles);
  void ReadIdFileSingle(int ifile, vector <Particle_t> &Particles);
  void ReadPosFileSingle(int ifile, vector <Particle_t> &Particles);
  void ReadVelFileSingle(int ifile, vector <Particle_t> &Particles);
  void CheckIdRange(vector <Particle_t> &Particles);
  void LoadExtraHeaderParams(JingHeader_t &header);
  void ReadParticleArray(float partarr[][3],long int np,int fileno);
public:
  JingReader_t(int snapshot_id): SnapshotId(snapshot_id)
  {
    NeedByteSwap=false;
    ProbeFiles();
  }
  void ReadHeader(JingHeader_t &header, const char filetype[]="pos", int ifile=0);
  void LoadSnapshot(vector <Particle_t> &Particles, Cosmology_t &Cosmology);
};


namespace JingGroup
{
  extern bool IsJingGroup(const string &GroupFileFormat);
  extern int ProbeGroupFileByteOrder(int snapshot_id);
  extern HBTInt LoadGroup(int snapshot_id, vector< Halo_t >& Halos);
};

#endif
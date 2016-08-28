#ifndef JING_IO_INCLUDED
#define JING_IO_INCLUDED

#include "../../halo.h"


#define OMEGA0 0.3  
#define OMEGAL0 0.7
#define PMass 1e-2
#define Redshift_INI 14.

struct JingHeader_t
{
   HBTInt Np;   //number of particles in the simulation
   HBTInt ips;  //time step number, or snapshot num
   HBTReal Redshift; // current redshift
   HBTReal OmegaM0;  //omega_m0
   HBTReal OmegaLambda0; //OmegaLambda0;   
   HBTReal Omegat;  // current omega_m, mass density
   HBTReal Lambdat; // current Omega_lambda, dark energy density
   HBTReal BoxSize;  // boxsize in unit of Mpc/h
   HBTReal xscale;
   HBTReal vscale;
   HBTReal Hz; //Hubble Param at ztp
   HBTReal vunit; //velocity unit to get physical peculiar velocity in km/s
   /*==extension needed for calculating binding energy:==*/
   HBTReal ScaleFactor;//current reduced scale factor 
   HBTReal mass[2];//Gas (mass[0]) and DM (mass[1]) particle masses, in units of 10^10Msun/h
};
class JingReader_t
{
  int SnapshotId;
  JingHeader_t Header;
  bool NeedByteSwap, FlagHasScale;
  int NumFilesPos, NumFilesVel, NumFilesId;
  
  string GetFileName(const char * filetype, int iFile=0);
  void CountFiles(const char *filetype, int &nfiles); 
  void ProbeFiles();
  void ReadId(vector <Particle_t> &Particles);
  void ReadVelocity(vector <Particle_t> &Particles);
  void ReadPosition(vector <Particle_t> &Particles);
  void ReadIdFileSingle(int ifile, vector <Particle_t> &Particles);
  void ReadPosFileSingle(int ifile, vector <Particle_t> &Particles);
  void ReadVelFileSingle(int ifile, vector <Particle_t> &Particles);
public:
  JingReader_t(int snapshot_id): SnapshotId(snapshot_id)
  {
    NeedByteSwap=false;
    ProbeFiles();
	FlagHasScale=(Nsnap<=SNAP_DIV_SCALE);//FIXME
  }
  void ReadHeader(JingHeader_t &header, const char filetype[]="pos", int ifile=0);
  void LoadSnapshot(vector <Particle_t> &Particles, Cosmology_t &Cosmology);
  void LoadGroup();
};

#endif
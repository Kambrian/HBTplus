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
   HBTReal ztp; // current redshift
   HBTReal Omega0;  //omega_m0
   HBTReal OmegaLambda; //OmegaLambda0;   
   HBTReal Omegat;  // current omega_m, mass density
   HBTReal Lambdat; // current Omega_lambda, dark energy density
   HBTReal rLbox;  // boxsize in unit of Mpc/h
   HBTReal xscale;
   HBTReal vscale;
   HBTReal Hz; //Hubble Param at ztp
   HBTReal vunit; //velocity unit to get physical peculiar velocity in km/s
   /*==extension needed for calculating binding energy:==*/
   HBTReal time;//current reduced scale factor 
   HBTReal mass[2];//Gas (mass[0]) and DM (mass[1]) particle masses, in units of 10^10Msun/h
};
class JingReader_t
{
  int SnapshotId;
  JingHeader_t Header;
  bool NeedByteSwap;
  int NumFilesPos, NumFilesVel, NumFilesId;
  
  string GetFileName(const char * filetype, int iFile=0);
  void CountFiles(const char *filetype, int &nfiles); 
  void ProbeFiles();
public:
  JingReader_t(int snapshot_id): SnapshotId(snapshot_id)
  {
    NeedByteSwap=false;
    ProbeFiles();
  }
  void ReadHeader(JingHeader_t &header);
  void LoadSnapshot(vector <Particle_t> &Particles, Cosmology_t &Cosmology);
  void LoadGroup();
};

#endif
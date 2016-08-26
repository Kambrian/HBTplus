typedef struct
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
   HBTInt Nsnap; //current snapnum
}IO_HEADER; //header of jing's data structure

typedef struct 
{
HBTInt Ngroups;
HBTInt Nids;
HBTInt *Len;
HBTInt *Offset;
HBTReal *HaloCen[3];//HaloCen[3][Ngroups], center of halos
HBTInt *PIDorIndex; //stores PID when loading FOF then changes to PIndex after loading particles
char *HaloMask; //HaloMask[NP_DM],HaloMask==1 means the particle does not belong to any sub,i.e, it's free.
char *HaloMaskSrc;
HBTInt *ID2Halo;//for index2halo, this is ID2halo[NP_DM];
}CATALOGUE;

struct ParticleData
{
//#ifndef PID_ORDERED
HBTInt *PID; //since Jing's snap data is stored according to PID, so this array is not necessary here
//#endif
HBTReal (*Pos)[3];
HBTReal (*Vel)[3];
HBTInt Nsnap; //this exists to help check whether the Pdat is loaded for the current snapshot
};	

extern HBTReal PMass;
extern IO_HEADER header;//this header filled during load_particle_data();
extern struct ParticleData Pdat;
extern HBTInt snaplist[MaxSnap];


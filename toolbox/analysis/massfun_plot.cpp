//outputs final massfunction data for plots
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "../../src/datatypes.h"
#include "../../src/config_parser.h"
#include "../../src/snapshot.h"
#include "../../src/halo.h"
#include "../../src/subhalo.h"
#include "../../src/mymath.h"
#include "../../src/linkedlist_parallel.h"

//#define NORM  //produce normalized massfunction (in terms Msub/Mhost rather than Msub)
#define RMIN 0
#define RMAX 1	 //statistics done in RMIN*rvi<r<RMAX*rvir
#define NBIN 20  //bin number for Msub
// #define AUTOBIN
//#define FIX_XBIN  //define this to use preset xmass bin
#define NFUN 5   //bin number for Mhost

#define MHOST M200Crit
#define RHOST R200CritComoving

float ParticleMass;
#define NBOUNDMIN 20

class MassFunc_t
{
private:
  vector <float> Mlist;  
  float XRange[2];
  bool HostInBin(float m)
  {
    return m>=Mbin[0]&&m<Mbin[1];
  }
public:  
  float Mbin[2]; //host mass range
  int Nhost;//number of hosts in this host mass bin
  double Mhost;//total host mass in this bin
  float XMassLow[NBIN], XMassMean[NBIN];//x coordinates,[Msub_lower_lim,Msub_av] for each Msub mass bin
  float MfunSpec[NBIN][2];//[dN/dMsub,sqrt(dN)/dMsub], massfun plus Poison error for each Msub bin
  float MfunSpecln[NBIN][2];//[dN/dlnMsub,...]
  float MfunCum[NBIN][2];//[N(>Msub_lower_lim),...]
  MassFunc_t():Mbin(), XRange(), Mlist(), Nhost(0), Mhost(0.), XMassLow{}, XMassMean{}, MfunSpec{}, MfunSpecln{}, MfunCum{}
  {}
  void build(float xrange[2], float mbin[2], const SubhaloSnapshot_t &subsnap, Linkedlist_t &ll)
  {
    Mbin[0]=mbin[0];
    Mbin[1]=mbin[1];
    XRange[0]=xrange[0];
    XRange[1]=xrange[1];
    mass_list(subsnap, ll);
    mass_count();
  }
  void mass_list(const SubhaloSnapshot_t &subsnap, Linkedlist_t &ll)
  {
    Mlist.reserve(subsnap.Subhalos.size());
    auto &subgroups=subsnap.MemberTable.SubGroups;
    for(HBTInt grpid=0;grpid<subgroups.size()-1;grpid++)
    {
      if(subgroups[grpid].size()==0) continue;
      HBTInt cenid=subgroups[grpid][0];
      auto &mhost=subsnap.Subhalos[cenid].MHOST;
      if(HostInBin(mhost))
      {	
	collect_submass(grpid,subsnap, ll);
	Nhost++;
	Mhost+=mhost;
      }
    }
  }
  void mass_count()
  {
    /* divide Mlist into nbin and count dN
     * xmass:size(nbin,2), [lower bin limits, center of mass bins]
     * MfunSpec: size(nbin,2), [dN/dMsub,error]
     * MfunSpecln: same as above for dN/dlnMsub
     * MfunCum: size(nbin,2), cumulative mass counts N(>Msub) 
     * */
    float xmin,xmax,dlnx;
    vector <double> mass(NBIN, 0.);
    vector <HBTInt> count(NBIN, 0);
    
    xmin=NBOUNDMIN*ParticleMass;
    #ifdef NORM
    float xmin_rel=NBOUNDMIN*ParticleMass/mfun->Mbin[0];
    if(xmin<xmin_rel) xmin=xmin_rel;
    #endif
    xmax=(*max_element(Mlist.begin(), Mlist.end()))*1.001;
    
    #ifdef FIX_XBIN
    xmin=XRange[0];
    xmax=XRange[1];
    #endif
    
    vector <float> x;
    logspace(xmin,xmax,NBIN+1, x);
    dlnx=logf(x[1]/x[0]);
    printf("%f,%f\n",xmin,xmax);
    
    for(auto &&m: Mlist)
    {
      int i=upper_bound(x.begin(), x.end(), m)-x.begin()-1;
      if(i>=0&&i<NBIN)
      {
	count[i]++;
	mass[i]+=m;
      }
    }
    
    HBTInt ncum=0;
    for(int i=NBIN-1;i>=0;i--)
    {
      ncum+=count[i];
      XMassLow[i]=x[i];
      XMassMean[i]=mass[i]/count[i];
      MfunSpec[i][0]=count[i]/(x[i+1]-x[i]);
      MfunSpec[i][1]=sqrt(count[i])/(x[i+1]-x[i]);
      MfunSpecln[i][0]=count[i]/dlnx;//another way: <m>*dN/dm=xmass[i,1]*MfunSpec;
      MfunSpecln[i][1]=sqrt(count[i])/dlnx;
      MfunCum[i][0]=ncum;
      MfunCum[i][1]=sqrt(ncum);
    }
  }
  
  
  void collect_submass(int grpid, const SubhaloSnapshot_t &subsnap, Linkedlist_t &ll)
  {
    if(subsnap.MemberTable.SubGroups[grpid].size()==0) return;
    HBTInt cenid=subsnap.MemberTable.SubGroups[grpid][0];
//     float rmin=RMIN*subsnap.Subhalos[cenid].RHOST;
    float rmax=RMAX*subsnap.Subhalos[cenid].RHOST;
#ifdef NORM
    float mhost=subsnap.Subhalos[cenid].MHOST;//TODO: try with bound mass??
#else
    float mhost=1.;
#endif
    HBTxyz &cenpos=subsnap.Subhalos[cenid].ComovingMostboundPosition;
    vector <HBTInt> sublist;
    ll.SearchSphere(rmax, cenpos, sublist, 8);
    Mlist.reserve(Mlist.size()+sublist.size());
    for(HBTInt i=0;i<sublist.size();i++)
    {
      if(sublist[i]!=cenid)
	Mlist.push_back(subsnap.Subhalos[i].Mbound/mhost);
    }
  }
  
  
  
  size_t save(hid_t file)
  {
    hsize_t ndim=1, dim_atom[]={1}, dims[2]={2,2};
    writeHDFmatrix(file, Mbin, "HostMassRange", 1, dims, H5T_NATIVE_FLOAT);
    writeHDFmatrix(file, Nhost, "NumberOfHosts", 1, dim_atom, H5T_NATIVE_INT);
    writeHDFmatrix(file, Mhost, "HostMassSum", 1, dim_atom, H5T_NATIVE_DOUBLE);
    dims[0]=NBIN;
    writeHDFmatrix(file, XMassLow, "SubMassLow", 1, dims, H5T_NATIVE_FLOAT);
    writeHDFmatrix(file, XMassLow, "SubMassMean", 1, dims, H5T_NATIVE_FLOAT);
    writeHDFmatrix(file, MfunSpec, "SpecificMF", 2, dims, H5T_NATIVE_FLOAT);
    writeHDFmatrix(file, MfunSpecln, "SpecificMFLn", 2, dims, H5T_NATIVE_FLOAT);
    writeHDFmatrix(file, MfunCum, "CumulativeMF", 2, dims, H5T_NATIVE_FLOAT);
  }
};

class SubhaloPos_t: public PositionData_t
{
  vector <Subhalo_t> &Subhalos;
public:
  SubhaloPos_t(vector <Subhalo_t> &subhalos):Subhalos(subhalos)
  {}
  const HBTxyz & operator [](HBTInt i) const
  {
    return Subhalos[i].ComovingMostboundPosition; 
  }
  size_t size() const
  {
    return Subhalos.size();
  }
};


void logspace(double xmin,double xmax,int N,vector <float> &x);
// float (*Mvir)[3],(*Rvir)[3];
// void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap);
int main(int argc,char **argv)
{
  int snapshot_start, snapshot_end, isnap;
  ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  isnap=HBTConfig.MaxSnapshotIndex;
  SubhaloSnapshot_t subsnap(isnap);
  ParticleMass=subsnap.Cosmology.ParticleMass;
  
  /* decide host mass bins */  
  vector <float> Mgrpbin={pow(10,2), pow(10, 2.5), pow(10,3), pow(10, 4), pow(10,5)};	  
  #ifdef AUTOBIN
  float Mmax=0.;
  auto &subgroups=subsnap.MemberTable.SubGroups;
  HBTInt Ngroups=subgroups.size()-1;
  for(HBTInt grpid=0;grpid<Ngroups;grpid++)
  {
    if(subgroups[grpid].size()==0) continue;
    HBTInt subid=subgroups[grpid][0];
    if(subsnap.Subhalos[subid].MHOST>Mmax)
      Mmax=subsnap.Subhalos[subid].MHOST;
  }
  Mmax=Mmax*1.01;
  float Mmin=1000*subsnap.Cosmology.ParticleMass; //if not def DM_ONLY, set it manually here
  logspace(Mmin,Mmax,NFUN+1,Mgrpbin);
  #endif
   
  float xrange[NFUN][2]={0.620373,50,
    0.620373,500,
    0.620373,5000};  //this only takes effect when FIX_XBIN is defined
    
    SubhaloPos_t SubPos(subsnap.Subhalos);
    LinkedlistPara_t ll(200, SubPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
    MassFunc_t mfun[NFUN];							
    #pragma omp parallel for
    for(int i=0;i<NFUN;i++)
      mfun.build(xrange[i], &Mgrpbin[i], subsnap, ll); 
    
    string outdir=HBTConfig.SubhaloPath+"/analysis/";
    mkdir(outdir.c_str(),0755);
    
    #ifdef NORM
    sprintf(buf,"%s/massfunN_%03d.%d",outputdir,Nsnap,VirType);
    #else
    sprintf(buf,"%s/massfun_%03d.%d",outputdir,Nsnap,VirType);
    #endif
    myfopen(fp,buf,"w");	
    {
      int ntmp;
      float rtmp;
      rtmp=1./header.time-1;
      fwrite(&rtmp,sizeof(float),1,fp);
      ntmp=NFUN;
      fwrite(&ntmp,sizeof(int),1,fp);
      ntmp=VirType;
      fwrite(&ntmp,sizeof(int),1,fp);
      rtmp=RMIN;
      fwrite(&rtmp,sizeof(int),1,fp);
      rtmp=RMAX;
      fwrite(&rtmp,sizeof(int),1,fp);
      //~ fwrite(mfun,sizeof(struct MassFunc),NFUN,fp);
      write_massfuncs(mfun,NFUN,fp);
      ntmp=NFUN;
      fwrite(&ntmp,sizeof(int),1,fp);
    }
    fclose(fp);
    myfree(Mvir);
    myfree(Rvir);
    
    return 0;
}

void logspace(double xmin,double xmax,int N, vector <float> &x)
{
  x.resize(N);
  int i;
  double dx;
  x[0]=xmin;x[N-1]=xmax;
  xmin=log(xmin);
  xmax=log(xmax);
  dx=exp((xmax-xmin)/(N-1));
  for(i=1;i<N-1;i++)
  {
    x[i]=x[i-1]*dx;
  }
}
/*
void load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
{
  char buf[1024];
  FILE *fp;
  int Nvir[3],i,j;
  sprintf(buf,"%s/profile/logbin/halo_size_%03d",SUBCAT_DIR,Nsnap);
  myfopen(fp,buf,"r");
  for(i=0;i<Ngroups;i++)
  {
    fseek(fp,14*4L,SEEK_CUR);
    fread(Nvir,sizeof(int),3,fp);
    for(j=0;j<3;j++)
      Mvir[i][j]=Nvir[j]*partmass;
    fread(Rvir+i,sizeof(float),3,fp);
    fseek(fp,4*4L,SEEK_CUR);
  }
  fclose(fp);
}*/
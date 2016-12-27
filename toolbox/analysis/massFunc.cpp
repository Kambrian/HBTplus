//outputs final massfunction data for plots
#include <cmath>
#include <iostream>
#include <string>

#include "../../src/datatypes.h"
#include "../../src/config_parser.h"
#include "../../src/snapshot.h"
#include "../../src/halo.h"
#include "../../src/subhalo.h"
#include "../../src/mymath.h"
#include "../../src/linkedlist_parallel.h"

#define NORM  //produce normalized massfunction (in terms Msub/Mhost rather than Msub)
// #define RMIN 0
#define RMAX 1	 //statistics done in RMIN*rvi<r<RMAX*rvir
#define NBIN 25  //bin number for Msub
#define FIX_HOSTBIN
// #define FIX_XBIN  //define this to use preset xmass bin
#define NFUN 7   //bin number for Mhost

#define EXTERN_VIR //whether to use external or internal virial
#define MHOST M200Crit
#define RHOST R200CritComoving
#define MSUB Mbound  //or LastMaxMass for unevolved MF

// #define EXCLUDE_EJECTED_HOST //exclude ejected halos from host list

// #define EXCLUDE_SOFTENING //exclude 10*SofteningHalo in the center

#define LIST_MAINSAT //list most-massive satellite of each host
//#define PARTICLEMASS 8.6e-2 //millimill
#define PARTICLEMASS 6.88e-4 //mill2
// #define PARTICLEMASS 0.000229431 #AqA5
float ParticleMass;

class MassFunc_t
{
private:
  vector <float> Mlist; 
#ifdef LIST_MAINSAT
  vector <float> MainSatList;//mass list of most-massive satellites
  vector <float> MainSatRList;
#endif
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
  float MfunSpecln[NBIN][2];//[dN/dlnMsub,...]
  float MfunCum[NBIN][2];//[N(>Msub_lower_lim),...]
  MassFunc_t():Mbin(), XRange(), Mlist(),  Nhost(0), Mhost(0.), XMassLow{}, XMassMean{}, MfunSpecln{}, MfunCum{}
#ifdef LIST_MAINSAT
, MainSatList(), MainSatRList()
#endif
  {}
  void build(float xrange[2], float mbin[2], const SubhaloSnapshot_t &subsnap, LinkedlistPara_t &ll);
  void mass_list(const SubhaloSnapshot_t &subsnap, LinkedlistPara_t &ll);
  void mass_count();  
  void collect_submass(int grpid, const SubhaloSnapshot_t &subsnap, LinkedlistPara_t &ll);
  void save(hid_t file);
};

class SubhaloPos_t: public PositionData_t
{
  vector <Subhalo_t> &Subhalos;
public:
  SubhaloPos_t(vector <Subhalo_t> &subhalos):Subhalos(subhalos)
  {}
  const HBTxyz & operator [](HBTInt i) const
  {    return Subhalos[i].ComovingMostBoundPosition; 
  }
  size_t size() const
  {    return Subhalos.size();
  }
};

void logspace(double xmin,double xmax,int N,vector <float> &x);
#ifdef EXTERN_VIR
struct HaloSize_t
{
  HBTInt HaloId;
  HBTxyz CenterComoving;
  float MVir, M200Crit, M200Mean, RVirComoving, R200CritComoving, R200MeanComoving;
  float RmaxComoving, VmaxPhysical;
};
vector <HaloSize_t> HaloSize;
void LoadHaloSize(vector<HaloSize_t> &HaloSize, int isnap);
#endif
int main(int argc,char **argv)
{
  if(argc!=3)
  {
    cerr<<"Usage: "<<endl;
    cerr<<" "<<argv[0]<<" [config_file] [snapshot_number]"<<endl;
    cerr<<"    If snapshot_number<0, then it's counted from final snapshot in reverse order"<<endl; 
    cerr<<"    (i.e., FinalSnapshot=-1,... FirstSnapshot=-N)"<<endl;
    return 1;
  }
  HBTConfig.ParseConfigFile(argv[1]);
  int isnap=atoi(argv[2]);
  if(isnap<0) isnap=HBTConfig.MaxSnapshotIndex+isnap+1;
  SubhaloSnapshot_t subsnap(isnap, SubReaderDepth_t::SubTable);
  ParticleMass=subsnap.Cosmology.ParticleMass;//if not def DM_ONLY, set it manually here
  if(ParticleMass==0)
  {
    ParticleMass=PARTICLEMASS;
    cerr<<"Error: failed to read particle mass in subhalo snapshot. Manually setting it to "
    <<PARTICLEMASS<<" according to macro."<<endl;
  }
  cout<<"data loaded\n";
#ifdef EXTERN_VIR
  LoadHaloSize(HaloSize,isnap);
  if(HaloSize.size()!=subsnap.MemberTable.SubGroups.size())
  {
    cerr<<HaloSize.size()<<"!="<<subsnap.MemberTable.SubGroups.size()<<endl;
    throw(runtime_error("number of halos mismatch\n"));
  }
#endif
  /* decide host mass bins */  
  float Mgrpbin[NFUN][2]={1e0, 1e1, 
    1e1, 1e2, 
    1e2, 1e3, 
    1e3, 1e4, 
    1e4, pow(10,4.5), 
    1e3, pow(10,3.5), 
    1e2, pow(10,2.5)};
  #ifndef FIX_HOSTBIN
  float Mmax=0.;
  auto &subgroups=subsnap.MemberTable.SubGroups;
  HBTInt Ngroups=subgroups.size();
#pragma omp parallel for reduction(max: Mmax)
  for(HBTInt grpid=0;grpid<Ngroups;grpid++)
  {
    if(subgroups[grpid].size()==0) continue;
    HBTInt cenid=subgroups[grpid][0];
#ifdef EXTERN_VIR
    auto &mhost=HaloSize[grpid].MHOST;
#else
    auto &mhost=subsnap.Subhalos[cenid].MHOST;
#endif
    if(mhost>Mmax)
      Mmax=mhost;
  }
  Mmax=Mmax*1.01;
  float Mmin=1000*ParticleMass; 
  vector <float> bin_edges(NFUN+1);
  logspace(Mmin,Mmax,NFUN+1,bin_edges);
  for(int i=0;i<NFUN;i++)
  {
    Mgrpbin[i][0]=bin_edges[i];
    Mgrpbin[i][1]=bin_edges[i+1];
  }
  #endif
  
#ifdef NORM
  float xrange[NFUN][2]={1e-6,1,
    1e-6,1,
    1e-6,1,
    1e-6,1,
    1e-6,1,
    1e-6,1,
    1e-6,1
  };  //this only takes effect when FIX_XBIN is defined
#else
  float xrange[NFUN][2]={0.620373,50,
    0.620373,500,
    0.620373,5000};
#endif
    
    SubhaloPos_t SubPos(subsnap.Subhalos);
    LinkedlistPara_t ll(200, &SubPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
    cout<<"linked list compiled\n";
    
    MassFunc_t mfun[NFUN];							
//     #pragma omp parallel for
    for(int i=0;i<NFUN;i++)
      mfun[i].build(xrange[i], Mgrpbin[i], subsnap, ll); 
    cout<<"mass func computed\n";
    
    string outdir=HBTConfig.SubhaloPath+"/analysis/";
    mkdir(outdir.c_str(),0755);
#ifdef EXCLUDE_EJECTED_HOST
    string suffix=".no_eject.hdf5";
#else
    string suffix=".hdf5";
#endif
#ifdef EXCLUDE_SOFTENING
    suffix=".softcut"+suffix;
#endif
#ifdef NORM
    string funcname="massFuncN";
#else
    string funcname="massFunc";
#endif
#define xstr(s) str(s)
#define str(s) #s
#ifdef EXTERN_VIR
    string filename=outdir+funcname+to_string(isnap)+"." xstr(MHOST) "." xstr(MSUB) +suffix;
#else
    string filename=outdir+funcname+to_string(isnap)+"." xstr(MHOST) "_bound." xstr(MSUB) +suffix;
#endif
#undef xstr
#undef str
    hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    for(int i=0;i<NFUN;i++)
    {
      string dsetname="/massFunc_"+to_string(i);
      hid_t grp=H5Gcreate2(file, dsetname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      mfun[i].save(grp);
      H5Gclose(grp);
    }
        
    return 0;
}

void MassFunc_t::build(float xrange[2], float mbin[2], const SubhaloSnapshot_t &subsnap, LinkedlistPara_t &ll)
{
  Mbin[0]=mbin[0];
  Mbin[1]=mbin[1];
  XRange[0]=xrange[0];
  XRange[1]=xrange[1];
  mass_list(subsnap, ll);
  mass_count();
}
void MassFunc_t::mass_list(const SubhaloSnapshot_t &subsnap, LinkedlistPara_t &ll)
{
  Mlist.reserve(subsnap.Subhalos.size());
  auto &subgroups=subsnap.MemberTable.SubGroups;
  for(HBTInt grpid=0;grpid<subgroups.size();grpid++)
  {
    if(subgroups[grpid].size()==0) continue;
    HBTInt cenid=subgroups[grpid][0];
#ifdef EXTERN_VIR
    auto &mhost=HaloSize[grpid].MHOST;
#else
    auto &mhost=subsnap.Subhalos[cenid].MHOST;
#endif
    if(HostInBin(mhost))
    {
#ifdef EXCLUDE_EJECTED_HOST
      if(subsnap.Subhalos[cenid].Mbound<subsnap.Subhalos[cenid].LastMaxMass*0.9) continue;
#endif
      collect_submass(grpid,subsnap, ll);
      Nhost++;
      Mhost+=mhost;
    }
  }
}
void MassFunc_t::mass_count()
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
  
  xmin=HBTConfig.MinNumPartOfSub*ParticleMass;
  #ifdef NORM
  xmin/=Mbin[0];
  #endif
  xmax=(*max_element(Mlist.begin(), Mlist.end()))*1.001;
  
  #ifdef FIX_XBIN
  xmin=XRange[0];
  xmax=XRange[1];
  #endif
  
  vector <float> x;
  logspace(xmin,xmax,NBIN+1, x);
  dlnx=logf(x[1]/x[0]);
  printf("%f,%f:  %f,%f\n",Mbin[0], Mbin[1], xmin,xmax);
  
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
    MfunSpecln[i][0]=count[i]/dlnx;//another way: <m>*dN/dm=xmass[i,1]*MfunSpec;
    MfunSpecln[i][1]=sqrt(count[i])/dlnx;
    MfunCum[i][0]=ncum;
    MfunCum[i][1]=sqrt(ncum);
  }
}


void MassFunc_t::collect_submass(int grpid, const SubhaloSnapshot_t &subsnap, LinkedlistPara_t &ll)
{
  if(subsnap.MemberTable.SubGroups[grpid].size()==0) return;
  HBTInt cenid=subsnap.MemberTable.SubGroups[grpid][0];
//   float rmin=RMIN*subsnap.Subhalos[cenid].RHOST;
#ifdef EXCLUDE_SOFTENING
  float rmin=10*HBTConfig.SofteningHalo;
#else
  float rmin=-1;
#endif
#ifdef EXTERN_VIR
  float rmax=RMAX*HaloSize[grpid].RHOST;
#else
  float rmax=RMAX*subsnap.Subhalos[cenid].RHOST;
#endif
  float mhost=1.;
  #ifdef NORM
  #ifdef EXTERN_VIR
  mhost=HaloSize[grpid].MHOST;
  #else
  mhost=subsnap.Subhalos[cenid].MHOST;//TODO: try with bound mass??
  #endif
  #endif
  auto &cenpos=subsnap.Subhalos[cenid].ComovingMostBoundPosition;
  vector <LocatedParticle_t> sublist;
  ll.SearchSphere(rmax, cenpos, sublist, 8, rmin);
  Mlist.reserve(Mlist.size()+sublist.size());
  float mmax=0., rmmax=0.;
  for(auto &&subid: sublist)
  {
    if(subid.id!=cenid)
    {
      float m=subsnap.Subhalos[subid.id].MSUB/mhost;
      Mlist.push_back(m);
      if(m>mmax)
      {
	mmax=m;
	rmmax=subid.d;
      }
    }
  }
#ifdef LIST_MAINSAT
  if(mmax>0) 
  {
    MainSatList.push_back(mmax);
    MainSatRList.push_back(rmmax);
  }
#endif
}



void MassFunc_t::save(hid_t file)
{
  hsize_t ndim=1, dim_atom[]={1}, dims[2]={2,2};
  writeHDFmatrix(file, Mbin, "HostMassRange", 1, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, &Nhost, "NumberOfHosts", 1, dim_atom, H5T_NATIVE_INT);
  writeHDFmatrix(file, &Mhost, "HostMassSum", 1, dim_atom, H5T_NATIVE_DOUBLE);
  dims[0]=NBIN;
  writeHDFmatrix(file, XMassLow, "SubMassLow", 1, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, XMassMean, "SubMassMean", 1, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, MfunSpecln, "SpecificMFLn", 2, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, MfunCum, "CumulativeMF", 2, dims, H5T_NATIVE_FLOAT);
#ifdef LIST_MAINSAT
  dims[0]=MainSatList.size();
  writeHDFmatrix(file, MainSatList.data(), "MainSatList", 1, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, MainSatRList.data(), "MainSatRList", 1, dims, H5T_NATIVE_FLOAT);
#endif
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

#ifdef EXTERN_VIR
void BuildHDFHaloSize(hid_t &H5T_dtype)
{
  H5T_dtype=H5Tcreate(H5T_COMPOUND, sizeof (HaloSize_t));
  hsize_t dims[2]={3,3};
  hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);

  #define InsertMember(x,t) H5Tinsert(H5T_dtype, #x, HOFFSET(HaloSize_t, x), t)//;cout<<#x<<": "<<HOFFSET(HaloSize_t, x)<<endl
  InsertMember(HaloId, H5T_HBTInt);
  InsertMember(CenterComoving, H5T_HBTxyz);
  InsertMember(R200CritComoving, H5T_NATIVE_FLOAT);
  InsertMember(R200MeanComoving, H5T_NATIVE_FLOAT);
  InsertMember(RVirComoving, H5T_NATIVE_FLOAT);
  InsertMember(M200Crit, H5T_NATIVE_FLOAT);
  InsertMember(M200Mean, H5T_NATIVE_FLOAT);
  InsertMember(MVir, H5T_NATIVE_FLOAT);
  InsertMember(RmaxComoving, H5T_NATIVE_FLOAT);
  InsertMember(VmaxPhysical, H5T_NATIVE_FLOAT);
  #undef InsertMember	

  H5Tclose(H5T_HBTxyz);
}
void LoadHaloSize(vector <HaloSize_t> &HaloSize, int isnap)
{
  string filename=HBTConfig.SubhaloPath+"/HaloSize/HaloSize_"+to_string(isnap)+".hdf5";
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t H5T_HaloSize;
  BuildHDFHaloSize(H5T_HaloSize);
  
  hid_t dset=H5Dopen2(file, "HostHalos", H5P_DEFAULT);
  hsize_t dims[1];
  GetDatasetDims(dset, dims);
  HaloSize.resize(dims[0]);
  if(dims[0])	H5Dread(dset, H5T_HaloSize, H5S_ALL, H5S_ALL, H5P_DEFAULT, HaloSize.data());
  H5Dclose(dset);
  
  H5Tclose(H5T_HaloSize);
  H5Fclose(file);
}
#endif
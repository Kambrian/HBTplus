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
#define RMIN 0
#define RMAX 1	 //statistics done in RMIN*rvi<r<RMAX*rvir
#define NBIN 25  //bin number for Msub
#define FIX_HOSTBIN
// #define FIX_XBIN  //define this to use preset xmass bin
#define NFUN 5   //bin number for Mhost

#define EXTERN_VIR //whether to use external or internal virial
#ifdef EXTERN_VIR
#define VirType "Crit200"
#else
#define MHOST M200Crit
#define RHOST R200CritComoving
#endif

#define MSUB Mbound  //or LastMaxMass for unevolved MF

//#define PARTICLEMASS 8.6e-2 //millimill
#define PARTICLEMASS 6.88e-4 //mill2
float ParticleMass;

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
  float MfunSpecln[NBIN][2];//[dN/dlnMsub,...]
  float MfunCum[NBIN][2];//[N(>Msub_lower_lim),...]
  MassFunc_t():Mbin(), XRange(), Mlist(), Nhost(0), Mhost(0.), XMassLow{}, XMassMean{}, MfunSpecln{}, MfunCum{}
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
  {
    return Subhalos[i].ComovingMostBoundPosition; 
  }
  size_t size() const
  {
    return Subhalos.size();
  }
};


void logspace(double xmin,double xmax,int N,vector <float> &x);
#ifdef EXTERN_VIR
vector <float>MVIR, RVIR;
template <class MyReal>
void LoadVirial(int SnapshotId, vector <MyReal> &Mvir, vector <MyReal> &Rvir, const string & virtype);
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
  LoadVirial(subsnap.GetSnapshotId(), MVIR, RVIR, VirType);
  if(MVIR.size()!=subsnap.MemberTable.SubGroups.size())
  {
    cerr<<MVIR.size()<<"!="<<subsnap.MemberTable.SubGroups.size()<<endl;
    throw(runtime_error("number of halos mismatch\n"));
  }
#endif
  /* decide host mass bins */  
  vector <float> Mgrpbin={pow(10, 0), pow(10,1), pow(10,2), pow(10,3), pow(10, 4), pow(10,4.5)};	  
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
    auto &mhost=MVIR[grpid];
#else
    auto &mhost=subsnap.Subhalos[cenid].MHOST;
#endif
    if(mhost>Mmax)
      Mmax=mhost;
  }
  Mmax=Mmax*1.01;
  float Mmin=1000*ParticleMass; 
  logspace(Mmin,Mmax,NFUN+1,Mgrpbin);
  #endif
  
#ifdef NORM
  float xrange[NFUN][2]={1e-6,1,
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
      mfun[i].build(xrange[i], &Mgrpbin[i], subsnap, ll); 
    cout<<"mass func computed\n";
    
    string outdir=HBTConfig.SubhaloPath+"/analysis/";
    mkdir(outdir.c_str(),0755);
#ifdef NORM
    string funcname="massFuncN";
#else
    string funcname="massFunc";
#endif
#define xstr(s) str(s)
#define str(s) #s
#ifdef EXTERN_VIR
    string filename=outdir+funcname+to_string(isnap)+"." VirType "." xstr(MSUB) ".hdf5";
#else
    string filename=outdir+funcname+to_string(isnap)+"." xstr(MHOST) "." xstr(MSUB) ".hdf5";
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
    auto &mhost=MVIR[grpid];
#else
    auto &mhost=subsnap.Subhalos[cenid].MHOST;
#endif
    if(HostInBin(mhost))
    {	
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
  //     float rmin=RMIN*subsnap.Subhalos[cenid].RHOST;
#ifdef EXTERN_VIR
  float rmax=RMAX*RVIR[grpid];
#else
  float rmax=RMAX*subsnap.Subhalos[cenid].RHOST;
#endif
  #ifdef NORM
  #ifdef EXTERN_VIR
  float mhost=MVIR[grpid];
  #else
  float mhost=subsnap.Subhalos[cenid].MHOST;//TODO: try with bound mass??
  #endif
  #else
  float mhost=1.;
  #endif
  auto &cenpos=subsnap.Subhalos[cenid].ComovingMostBoundPosition;
  vector <HBTInt> sublist;
  ll.SearchSphere(rmax, cenpos, sublist);
  Mlist.reserve(Mlist.size()+sublist.size());
  for(auto &&subid: sublist)
  {
    if(subid!=cenid)
      Mlist.push_back(subsnap.Subhalos[subid].MSUB/mhost);
  }
}



void MassFunc_t::save(hid_t file)
{
  hsize_t ndim=1, dim_atom[]={1}, dims[2]={2,2};
  writeHDFmatrix(file, Mbin, "HostMassRange", 1, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, &Nhost, "NumberOfHosts", 1, dim_atom, H5T_NATIVE_INT);
  writeHDFmatrix(file, &Mhost, "HostMassSum", 1, dim_atom, H5T_NATIVE_DOUBLE);
  dims[0]=NBIN;
  writeHDFmatrix(file, XMassLow, "SubMassLow", 1, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, XMassLow, "SubMassMean", 1, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, MfunSpecln, "SpecificMFLn", 2, dims, H5T_NATIVE_FLOAT);
  writeHDFmatrix(file, MfunCum, "CumulativeMF", 2, dims, H5T_NATIVE_FLOAT);
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
 v oi*d load_halo_virial_size(float Mvir[][3],float Rvir[][3],float partmass,int Ngroups,int Nsnap)
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
#ifdef EXTERN_VIR
#include "../../src/io/gadget_group_io.h"
#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
template <class MyReal>
void LoadVirial(int SnapshotId, vector <MyReal> &Mvir, vector <MyReal> &Rvir, const string & virtype)
{
//   typedef float MyReal;
  
  int itype;
  if(virtype=="Mean200")
    itype=0;
  else if(virtype=="Crit200")
    itype=1;
  else if(virtype=="TopHat")
    itype=2;
  else
    throw(runtime_error("unknow virtype"+virtype));
	
  FILE *fd;
  char filename[1024];
  vector <int> Len;
  vector <HBTInt> Offset;
  int Ngroups, TotNgroups, Nids, NFiles;
  HBTInt NumberOfHaloes;
  long long TotNids;
  bool NeedByteSwap;
  bool IsGroupV3=("gadget3_int"==HBTConfig.GroupFileFormat||"gadget3_long"==HBTConfig.GroupFileFormat);
  
  string filename_format;
  bool IsSubFile;
  int FileCounts;
  GadgetGroup::GetFileNameFormat(SnapshotId, filename_format, FileCounts, IsSubFile, NeedByteSwap);
  assert(IsSubFile);
  
  long long Nload=0;
  for(int iFile=0;iFile<FileCounts;iFile++)
  {
	if(FileCounts>1)
	  sprintf(filename, filename_format.c_str(), "tab", iFile);
	else
	  sprintf(filename, filename_format.c_str(), "tab");
		
	myfopen(fd,filename,"r");
	myfread(&Ngroups, sizeof(Ngroups), 1, fd);
	if(IsGroupV3)
	{
	  myfread(&TotNgroups, sizeof(TotNgroups), 1, fd);
	  myfread(&Nids, sizeof(Nids), 1, fd);
	  myfread(&TotNids,sizeof(TotNids),1,fd);
	}
	else
	{
	  myfread(&Nids, sizeof(Nids), 1, fd);
	  myfread(&TotNgroups, sizeof(TotNgroups), 1, fd);
	}
	myfread(&NFiles, sizeof(NFiles), 1, fd);
	int Nsub,TotNsub;
	myfread(&Nsub,sizeof(int),1,fd);
	myfread(&TotNsub,sizeof(int),1,fd);
	if(FileCounts!=NFiles)
	{
	  cout<<"File count mismatch for file "<<filename<<": expect "<<FileCounts<<", got "<<NFiles<<endl;
	  exit(1);
	}
	
	if(0==iFile)
	{
	  NumberOfHaloes=TotNgroups; 	  
	  Mvir.resize(TotNgroups);
	  Rvir.resize(TotNgroups);
	}
	
	if(IsGroupV3)
	fseek(fd, (sizeof(int)*2+sizeof(MyReal)*4)*Ngroups, SEEK_CUR);//skip len,offset,mass,pos
	else
	  fseek(fd, sizeof(int)*(2*Ngroups+3*Nsub), SEEK_CUR); 
	fseek(fd, sizeof(MyReal)*Ngroups*2*itype, SEEK_CUR);
	myfread(Mvir.data()+Nload, sizeof(MyReal), Ngroups, fd);
	myfread(Rvir.data()+Nload, sizeof(MyReal), Ngroups, fd);
	if(feof(fd))
	{
	  fprintf(stderr,"error:End-of-File in %s\n",filename);
	  fflush(stderr);exit(1);  
	}
	Nload+=Ngroups;
	fclose(fd);
  }
  if(Nload!=NumberOfHaloes)
  {
	cerr<<"error:Num groups loaded not match: "<<Nload<<','<<NumberOfHaloes<<"\n";
	exit(1);
  }
  cout<<"Halo size "<<virtype<<"loaded at snapshot "<<SnapshotId<<", TotNgroups="<<TotNgroups<<" NFiles="<<NFiles<<endl;
  
}
#endif
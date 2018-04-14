//to compute the virial sizes, rmax and vmax of the host halos, using the ComovingMostBoundPosition of central subhalos as reference center
#include <cmath>
#include <iostream>
#include <string>

#include "../src/datatypes.h"
#include "../src/config_parser.h"
#include "../src/snapshot.h"
#include "../src/halo.h"
#include "../src/subhalo.h"
#include "../src/mymath.h"
#include "../src/linkedlist.h"

#define RMAX 2.

// #define NBIN 10  //define this to also compute a count profile; the first bin is [0, 1e-2)*Rvir and the remaining bins are logspaced between [1e-2,1)*Rvir
#ifdef NBIN
float Dlnx=(0-(-2.))/(NBIN-1); //1e-2 to 1, logspace
#endif

#define NDIV 256

struct HaloSize_t
{
  HBTInt HaloId;
  HBTxyz CenterComoving;
  float MVir, M200Crit, M200Mean, RVirComoving, R200CritComoving, R200MeanComoving;
  float RmaxComoving, VmaxPhysical;
#ifdef NBIN
  HBTInt Profile[NBIN];
#endif
  HaloSize_t(){};
  void fill_zeros()
  {
    CenterComoving[0]=CenterComoving[1]=CenterComoving[2]=0.;
    MVir=0.;
    M200Crit=0.;
    M200Mean=0.;
    RVirComoving=0.;
    R200CritComoving=0.;
    R200MeanComoving=0.;
    RmaxComoving=0.;
    VmaxPhysical=0.;
#ifdef NBIN
    fill(begin(Profile), end(Profile), 0);
#endif
  }
  void Compute(HBTxyz &cen, float rmax, HBTInt nguess, Linkedlist_t &ll, const ParticleSnapshot_t &partsnap);
};
void BuildHDFHaloSize(hid_t &H5T_dtypeInMem, hid_t &H5T_dtypeInDisk);
inline bool CompProfRadiusVal(const RadMassVel_t &a, const float r)
{
  return a.r<r;
}
inline bool CompProfRadius(const RadMassVel_t &a, const RadMassVel_t &b)
{
  return a.r<b.r;
}
inline bool CompProfVel(const RadMassVel_t &a, const RadMassVel_t &b)
{
  return a.v<b.v;
}
inline float ComovingMean200Radius(float M200b, float OmegaM0)
{
  return pow(PhysicalConst::G*M200b/100./OmegaM0/PhysicalConst::H0/PhysicalConst::H0, 1./3); //physical R200b at z=0, is comoving at all redshift.
}
HBTReal virialF_tophat, virialF_b200, virialF_c200;
HBTReal VelocityUnit;
void save(vector <HaloSize_t> &HaloSize, int isnap);
int main(int argc, char **argv)
{
  /*
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
  */
  int snap_start, snap_end;
  ParseHBTParams(argc, argv, HBTConfig, snap_start, snap_end);
  for(int isnap=snap_start;isnap<=snap_end;isnap++)
  {
    Timer_t timer;
    timer.Tick();
    HaloSnapshot_t halosnap(isnap);
    SubhaloSnapshot_t subsnap(isnap, SubReaderDepth_t::SubTable);;
    ParticleSnapshot_t partsnap(isnap, false);
    auto &Cosmology=partsnap.Cosmology;
    Cosmology.HaloVirialFactors(virialF_tophat, virialF_b200, virialF_c200);
    VelocityUnit=PhysicalConst::G/partsnap.Cosmology.ScaleFactor;
    timer.Tick();cout<<"load: "<<timer.GetSeconds(1)<<" seconds\n";
  
    SnapshotPos_t PartPos(partsnap);
    Linkedlist_t ll(NDIV, &PartPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
    cout<<"linked list compiled\n";
    timer.Tick();cout<<"link: "<<timer.GetSeconds(2)<<" seconds\n";
    vector <HaloSize_t> HaloSize(halosnap.size());
    #pragma omp parallel for schedule(dynamic, 1)
    for(HBTInt grpid=0;grpid<halosnap.size();grpid++)
    {
      HaloSize[grpid].HaloId=grpid;//need to adjust this in MPI version..
      auto &subgroup=subsnap.MemberTable.SubGroups[grpid];
      if(subgroup.size()==0)
      {
	HaloSize[grpid].fill_zeros();
	continue;
      }
      HBTInt np=halosnap.Halos[grpid].size();
      float rmax=ComovingMean200Radius(np*Cosmology.ParticleMass, Cosmology.OmegaM0)*RMAX;//use b200 as a ref
      HaloSize[grpid].Compute(subsnap.Subhalos[subsnap.MemberTable.SubGroups[grpid][0]].ComovingMostBoundPosition, rmax, np, ll, partsnap);   
    }
    timer.Tick();cout<<"compute: "<<timer.GetSeconds(3)<<" seconds\n";
    save(HaloSize, isnap);
  }
  
  return 0;
}
class ProfCollector_t: public ParticleCollector_t
{
  const ParticleSnapshot_t &PartSnap;
public:
  vector <RadMassVel_t> Prof;
  ProfCollector_t(const ParticleSnapshot_t &partsnap, HBTInt n_reserve=0):PartSnap(partsnap), Prof()
  {
    Prof.reserve(n_reserve);
  }
  void Collect(HBTInt index, HBTReal d2)
  {
    Prof.emplace_back(sqrt(d2), PartSnap.GetMass(index));
  }
};
void HaloSize_t::Compute(HBTxyz &cen, float rmax, HBTInt nguess, Linkedlist_t &ll, const ParticleSnapshot_t &partsnap)
{
    copyHBTxyz(CenterComoving, cen);
    ProfCollector_t collector(partsnap, nguess);
    ll.SearchSphere(rmax, cen, collector);
    HBTInt np=collector.Prof.size();
    auto &prof=collector.Prof;
    sort(prof.begin(), prof.end(), CompProfRadius);
    double m_cum=0.;
    for(auto && p: prof)  p.m=(m_cum+=p.m);
    for(HBTInt i=0;i<np;i++)
    {
	    if(prof[i].r<HBTConfig.SofteningHalo) prof[i].r=HBTConfig.SofteningHalo; //resolution
	    prof[i].v=prof[i].m/prof[i].r;//v**2
    }
   
    partsnap.Cosmology.SphericalOverdensitySize(MVir, RVirComoving, virialF_tophat, prof);
    partsnap.Cosmology.SphericalOverdensitySize(M200Crit, R200CritComoving, virialF_c200, prof);
    partsnap.Cosmology.SphericalOverdensitySize(M200Mean, R200MeanComoving, virialF_b200, prof);
    
    auto it_rv=lower_bound(prof.begin(), prof.end(), R200CritComoving, CompProfRadiusVal);
    auto maxprof=max_element(prof.begin(), it_rv, CompProfVel);//restrict Rmax<R200C
    RmaxComoving=maxprof->r;
    VmaxPhysical=sqrt(maxprof->v*VelocityUnit);
    
#ifdef NBIN
    fill(begin(Profile), end(Profile), 0);
    float r0=1e-2*RVirComoving;
    for(auto &&p: prof)
    {
      if(p.r<RVirComoving)
      {
	float logr=log10f(p.r/r0);
	int ibin=ceilf(log10f(p.r/r0)/Dlnx); //1 to NBIN-1 from r0 to rvir
	if(ibin<0) 
	  ibin=0;
	else if(ibin>=NBIN) 
	  ibin=NBIN-1;
	Profile[ibin]++;
      }
    }
#endif
}

void save(vector <HaloSize_t> &HaloSize, int isnap)
{
  string filepath=HBTConfig.SubhaloPath+"/HaloSize";
  stringstream formatter;
  formatter<<filepath<<"/HaloSize_"<<setw(3)<<setfill('0')<<isnap<<".hdf5";
  string filename=formatter.str();

  mkdir(filepath.c_str(), 0755);
  hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t ndim=1, dim_atom[]={1}, dims[2]={HaloSize.size(),2};
  hid_t H5T_HaloSizeInMem, H5T_HaloSizeInDisk;
  BuildHDFHaloSize(H5T_HaloSizeInMem, H5T_HaloSizeInDisk);
  writeHDFmatrix(file, HaloSize.data(), "HostHalos", 1, dims, H5T_HaloSizeInMem, H5T_HaloSizeInDisk);
  H5Tclose(H5T_HaloSizeInDisk);
  H5Tclose(H5T_HaloSizeInMem);
  H5Fclose(file);
}

void BuildHDFHaloSize(hid_t &H5T_dtypeInMem, hid_t &H5T_dtypeInDisk)
{
  H5T_dtypeInMem=H5Tcreate(H5T_COMPOUND, sizeof (HaloSize_t));
  hsize_t dims[2]={3,3};
  hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);
//   hid_t H5T_FloatVec3=H5Tarray_create2(H5T_NATIVE_FLOAT, 1, dims);

  #define InsertMember(x,t) H5Tinsert(H5T_dtypeInMem, #x, HOFFSET(HaloSize_t, x), t)//;cout<<#x<<": "<<HOFFSET(HaloSize_t, x)<<endl
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
#ifdef NBIN
  dims[0]=NBIN;
  hid_t H5T_HBTIntVec=H5Tarray_create2(H5T_HBTInt, 1, dims);
  InsertMember(Profile, H5T_HBTIntVec);
  H5Tclose(H5T_HBTIntVec);
#endif
  #undef InsertMember	

  H5T_dtypeInDisk=H5Tcopy(H5T_dtypeInMem);
  H5Tpack(H5T_dtypeInDisk); //clear fields not added to save disk space

//   H5Tclose(H5T_FloatVec3);
  H5Tclose(H5T_HBTxyz);
}

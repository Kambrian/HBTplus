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
#include "../src/linkedlist_parallel.h"

#define RMAX 2.

struct HaloSize_t
{
  HBTInt HaloId;
  HBTxyz CenterComoving;
  float MVir, M200Crit, M200Mean, RVirComoving, R200CritComoving, R200MeanComoving;
  float RmaxComoving, VmaxPhysical;
  void Compute(HBTxyz &cen, float rmax, HBTInt nguess, LinkedlistPara_t &ll, ParticleSnapshot_t &partsnap);
};
void BuildHDFHaloSize(hid_t &H5T_dtypeInMem, hid_t &H5T_dtypeInDisk);

inline bool CompProfRadius(const RadVelMass_t &a, const RadVelMass_t &b)
{
  return a.r<b.r;
}
inline bool CompProfVel(const RadVelMass_t &a, const RadVelMass_t &b)
{
  return a.v<b.v;
}
inline float ComovingMean200Radius(float M200b, float OmegaM0)
{
  return pow(PhysicalConst::G*M200b/100./OmegaM0/PhysicalConst::H0/PhysicalConst::H0, 1./3); //physical R200b at z=0, is comoving at all redshift.
}
HBTReal virialF_tophat, virialF_b200, virialF_c200;
HBTReal VelocityUnit;
void save(vector <HaloSize_t> HaloSize, int isnap);
int main(int argc, char **argv)
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
  
  HaloSnapshot_t halosnap(isnap);
  SubhaloSnapshot_t subsnap(isnap, SubReaderDepth_t::SubTable);;
  ParticleSnapshot_t partsnap(isnap);
  auto &Cosmology=partsnap.Cosmology;
  Cosmology.HaloVirialFactors(virialF_tophat, virialF_b200, virialF_c200);
  VelocityUnit=PhysicalConst::G/partsnap.Cosmology.ScaleFactor;
 
  SnapshotPos_t PartPos(partsnap);
  LinkedlistPara_t ll(256, &PartPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
  cout<<"linked list compiled\n";
  
  vector <HaloSize_t> HaloSize(halosnap.size());
  #pragma omp parallel for
  for(HBTInt grpid=0;grpid<halosnap.size();grpid++)
  {
    HaloSize[grpid].HaloId=grpid;//need to adjust this in MPI version..
    auto &subgroup=subsnap.MemberTable.SubGroups[grpid];
    if(subgroup.size()==0)
      continue;
    HBTInt np=halosnap.Halos[grpid].size();
    float rmax=ComovingMean200Radius(np*Cosmology.ParticleMass, Cosmology.OmegaM0)*RMAX;//use b200 as a ref
    HaloSize[grpid].Compute(subsnap.Subhalos[subsnap.MemberTable.SubGroups[grpid][0]].ComovingMostBoundPosition, rmax, np, ll, partsnap);   
  }
  save(HaloSize, isnap);
}


void HaloSize_t::Compute(HBTxyz &cen, float rmax, HBTInt nguess, LinkedlistPara_t &ll, ParticleSnapshot_t &partsnap)
{
    copyHBTxyz(CenterComoving, cen);
    vector <LocatedParticle_t> founds;
    founds.reserve(nguess);
    ll.SearchSphereSerial(rmax, cen, founds);
    HBTInt np=founds.size();
    vector <RadVelMass_t> prof(np);
    for(HBTInt i=0;i<np;i++)
    {
	prof[i].r=sqrt(founds[i].d2);
	prof[i].m=partsnap.GetMass(founds[i].index);
    }
    sort(prof.begin(), prof.end(), CompProfRadius);
    double m_cum=0.;
    for(auto && p: prof)  p.m=(m_cum+=p.m);
    for(HBTInt i=0;i<np;i++)
    {
	    if(prof[i].r<HBTConfig.SofteningHalo) prof[i].r=HBTConfig.SofteningHalo; //resolution
	    prof[i].v=prof[i].m/prof[i].r;//v**2
    }
    auto maxprof=max_element(prof.begin(), prof.end(), CompProfVel);
    RmaxComoving=maxprof->r;
    VmaxPhysical=sqrt(maxprof->v*VelocityUnit);
   
    partsnap.Cosmology.SphericalOverdensitySize(MVir, RVirComoving, virialF_tophat, prof);
    partsnap.Cosmology.SphericalOverdensitySize(M200Crit, R200CritComoving, virialF_c200, prof);
    partsnap.Cosmology.SphericalOverdensitySize(M200Mean, R200MeanComoving, virialF_b200, prof);
}

void save(vector <HaloSize_t> HaloSize, int isnap)
{
  string filename=HBTConfig.SubhaloPath+"/HaloSize";
  mkdir(filename.c_str(), 0755);
  filename=filename+"/HaloSize_"+to_string(isnap)+".hdf5";
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
  #undef InsertMember	

  H5T_dtypeInDisk=H5Tcopy(H5T_dtypeInMem);
  H5Tpack(H5T_dtypeInDisk); //clear fields not added to save disk space

//   H5Tclose(H5T_FloatVec3);
  H5Tclose(H5T_HBTxyz);
}
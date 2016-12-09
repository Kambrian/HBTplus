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

#define RMAX 2.

class ParticlePos_t: public PositionData_t
{
  const ParticleSnapshot_t &Snap;
public:
  SubhaloPos_t(const ParticleSnapshot_t &snap):Snap(snap)
  {}
  const HBTxyz & operator [](HBTInt i) const
  {
    return Snap.GetComovingPosition(i);
  }
  size_t size() const
  {
    return Snap.size();
  }
};

struct HaloSize_t
{
  HBTxyz CenterComoving;
  float MVir, M200Crit, M200Mean, RVirComoving, R200CritComoving, R200MeanComoving;
  float RmaxComoving, VmaxPhysical;
  void Compute(HBTxyz &cen, HBTxyz rmax, HBTInt nguess, LinkedlistPara_t &ll, ParticleSnapshot_t &partsnap);
};

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
 
  ParticlePos_t PartPos(partsnap);
  LinkedlistPara_t ll(128, &PartPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
  cout<<"linked list compiled\n";
  
  vector <HaloSize_t> HaloSize(halosnap.size());
  for(HBTInt grpid=0;grpid<halosnap.size();grpid++)
  {
    auto &subgroup=subsnap.MemberTable.SubGroups[grpid];
    if(subgroup.size()==0)
      continue;
    HBTInt np=halosnap.Halos[grpid].size();
    float rmax=ComovingMean200Radius(np*Cosmology.ParticleMass, Cosmology.OmegaM0)*RMAX;//use b200 as a ref
    HaloSize[grpid].Compute(subsnap.Subhalos[subsnap.MemberTable.SubGroups[grpid][0]].ComovingMostBoundPosition, rmax, np, ll, partsnap);   
  }
  save(HaloSize, isnap);
}


void HaloSize_t::Compute(HBTxyz &cen, HBTxyz rmax, HBTInt nguess, LinkedlistPara_t &ll, ParticleSnapshot_t &partsnap)
{
    copyHBTxyz(CenterComoving, cen);
    vector <LocatedParticle_t> founds;
    ll.SearchSphere(rmax, cen, founds, nguess);
    HBTInt np=founds.size();
    vector <RadVelMass_t> prof(np);
    for(HBTInt i=0;i<np;i++)
    {
	prof[i].r=founds[i].d;
	prof[i].m=partsnap.GetMass(founds[i].id);
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
  string filename=HBTConfig.SubhaloPath+"/halosize";
  mkdir(filename.c_str(), 0755);
  filename=filename+"/HaloSize_"+to_string(isnap)+".hdf5";
  hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
}
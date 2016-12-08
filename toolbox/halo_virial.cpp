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
 
  ParticlePos_t PartPos(partsnap);
  LinkedlistPara_t ll(128, &partsnap, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
  cout<<"linked list compiled\n";
    
  for(HBTInt grpid=0;grpid<halosnap.size();grpid++)
  {
    auto &subgroup=subsnap.MemberTable.SubGroups[grpid];
    if(subgroup.size()==0)
    {
      init....
      continue;
    }
    rest=estimate_grp_size();//use b200 as a ref
    vector <LocatedParticle_t> founds;
    ll.SearchSphere(rest*RMAX, subsnap.Subhalos[subgroup[0]], founds, halosnap.Halos[grpid].size());
    vector <HBTReal> d(founds.size());
    for(HBTInt i=0;i<founds.size();i++) d[i]=founds[i].d;
    sort(d.begin(),d.end());
    partsnap.SphericalOverdensitySize(....);
  }
}
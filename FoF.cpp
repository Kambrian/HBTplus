using namespace std;
#include <iostream>
#include <string>
#include <cstdlib>
#include <omp.h>

#include "src/datatypes.h"
#include "src/config_parser.h"
#include "src/snapshot.h"
#include "src/io/hbt_group_io.h"
#include "src/linkedlist.h"
#include "src/fof_builder.h"

#define TREE_FOF 1
#define LL_FOF 2
/*Notes on fof performance:
 * Tree search scales better than LinkedList as clustering increases.
 * LinkedList search is faster when the cell size is smaller, maybe saturates when each cell contains ~1 particle. 
 * For weakly clustered data (high redshift), linkedlist can be faster than tree search; but this is useless, since tree search is already fast enough at high redshift as well.
 * In highly clustered region, the spatial complexity (memory consumption) can be quite high when trying to reach 1 particle cell size. While the tree code spatial complexity is acceptable (~2times particles).
 * This means tree search is perferred over LL.
*/
#ifndef FOF_METHOD 
#define FOF_METHOD TREE_FOF
#endif

#ifndef DM_ONLY
  #error "fof currently only works on a single particle type"
#endif

typedef vector <HBTInt> halo_t;
void build_group_catalogue(HBTReal linklength, Snapshot_t &snapshot, vector <halo_t> & halos);

int main(int argc,char **argv)
{
#ifdef _OPENMP
omp_set_nested(0);
#endif
int snapshot_start, snapshot_end, isnap;
ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
mkdir(HBTConfig.HaloPath.c_str(), 0755);
HBTConfig.DumpParameters(HBTConfig.HaloPath);

for(isnap=snapshot_start; isnap<=snapshot_end; ++isnap)
{
  ParticleSnapshot_t snapshot(isnap, false);
  double mean_density=HBTConfig.FoFDarkMatterMassFraction*snapshot.Cosmology.OmegaM0*3.*PhysicalConst::H0*PhysicalConst::H0/8./M_PI/PhysicalConst::G;
  double ParticleMass=snapshot.Cosmology.ParticleMass;
  HBTReal r=HBTConfig.FoFLinkParam*pow(ParticleMass/mean_density,1.0/3.0);
  Timer_t timer;
  timer.Tick();
  vector <halo_t> halos;
  build_group_catalogue(r, snapshot, halos);
  timer.Tick();
  cout<<timer.GetSeconds(1)<<" seconds for snapshot"<<isnap<<endl;

  // load_particle_data_bypart(Nsnap,SNAPSHOT_DIR,FLAG_LOAD_ID);//only load id
  /*for(auto &&h: halos)//has been converted in build_group_catalogue
  {
    for(auto &&p: h)
      p=snapshot.GetParticleId(p);
  }*/
  HBTGroupIO::SaveHDFGroups(isnap, snapshot.GetSnapshotId(), halos);
}

return 0;
}

struct GrpIDCompor_t
{
  vector<HBTInt> &TagLen;
  GrpIDCompor_t(vector <HBTInt> &grplen): TagLen(grplen)
  {
  }
  bool operator()(HBTInt i, HBTInt j)
  {
    return TagLen[i]>TagLen[j];
  }
};

void build_group_catalogue(HBTReal linklength, Snapshot_t &snapshot, vector <halo_t> &halos)
{
#if FOF_METHOD==TREE_FOF
FoFBuilder_t builder(linklength, snapshot);
builder.Link();
auto &TagLen=builder.GrpLen;
auto &ParticleTags=builder.GrpTags;
#endif
#if FOF_METHOD==LL_FOF
vector <HBTInt> TagLen, ParticleTags;
int ndiv=256;
LinkedlistLinkGroup(linklength, snapshot, TagLen, ParticleTags, ndiv);
#endif
//sort tags according to taglen
vector <HBTInt> SortedTags(TagLen.size());
for(HBTInt i=0;i<SortedTags.size();i++)
  SortedTags[i]=i;
GrpIDCompor_t compor(TagLen);
sort(SortedTags.begin(), SortedTags.end(), compor);
//revert tag to index
vector <HBTInt> tag2id(SortedTags.size());
for(HBTInt i=0;i<SortedTags.size();i++)
  tag2id[SortedTags[i]]=i;

halos.resize(TagLen.size());
HBTInt ngrp, grplen_min=HBTConfig.FoFSizeMin;
for(ngrp=0;ngrp<SortedTags.size();ngrp++)
{
  auto n=TagLen[SortedTags[ngrp]];
  if(n<grplen_min) break;
  halos[ngrp].reserve(n);
}

for(HBTInt i=0;i<ParticleTags.size();i++)
{
  auto id=tag2id[ParticleTags[i]];
  if(id<ngrp)
    halos[id].push_back(snapshot.GetMemberId(i));//these are particle ids
}
halos.resize(ngrp);

cout<<halos.size()<<" groups above mass limit: ";
for(int i=0;i<3;i++)
{
  if(i<halos.size())
    cout<<halos[i].size()<<", ";
}
cout<<"...\n";
}
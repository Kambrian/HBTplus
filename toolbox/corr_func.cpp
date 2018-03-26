/*to compute the average count profile around randomly selected particles
 */

#include <cmath>
#include <iostream>
#include <string>
#include <random>

#include "../src/datatypes.h"
#include "../src/config_parser.h"
#include "../src/snapshot.h"
#include "../src/halo.h"
#include "../src/subhalo.h"
#include "../src/mymath.h"
#include "../src/linkedlist_parallel.h"
#include "../src/geometric_tree.h"

#define SAMPLESIZE 100000
#define RMIN 1e-2 //only for deciding binwidth; not used to set the innermost bin edge.
#define RMAX 2e1
#define NBIN 20
// #define USE_LL //algorithm: whether to use linkedlist or geotree for spatial search.
#define NLOOP 1

class SnapshotOffset_t: public Snapshot_t
{
public:
  HBTInt Offset;
  HBTInt N;
  Snapshot_t & Snapshot;
  SnapshotOffset_t(HBTInt offset, Snapshot_t & fullsnapshot): Offset(offset), Snapshot(fullsnapshot), Snapshot_t(fullsnapshot)
  {
    N=Snapshot.size()-Offset;
  }
  HBTInt ReSize(HBTInt n)
  {
    N=n;
  }
  HBTInt size() const
  {
	return N;
  }
  HBTInt GetMemberId(const HBTInt i) const
  {
	return i+Offset;
  }
  HBTReal GetMass(const HBTInt i) const
  {
	return Snapshot.GetMass(GetMemberId(i));
  }
  const HBTxyz & GetPhysicalVelocity(const HBTInt i) const
  {
	return Snapshot.GetPhysicalVelocity(GetMemberId(i));
  }
  const HBTxyz & GetComovingPosition(const HBTInt i) const
  {
	return Snapshot.GetComovingPosition(GetMemberId(i));
  }
};

struct HaloSize_t
{
  HBTInt n[NBIN];
  HaloSize_t(){};
  void Compute(const HBTxyz &cen, LinkedlistPara_t &ll);
  void Compute(const HBTxyz &cen, GeoTree_t &tree);
};

float DlnX=logf(RMAX/RMIN)/(NBIN-1)*2; //spacing in ln-space for r**2
void save(float prof[], int isnap);
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
  
  SubhaloSnapshot_t subsnap(isnap, SubReaderDepth_t::SubTable);;
  ParticleSnapshot_t partsnap(isnap);
  auto &Cosmology=partsnap.Cosmology;
    
  vector <HaloSize_t> Profiles(SAMPLESIZE);
#pragma omp parallel for
  for(int i=0;i<SAMPLESIZE;i++)
    fill(begin(Profiles[i].n), end(Profiles[i].n), 0);
    
  const int nloop=NLOOP; //do it in blocks to save memory
  HBTInt nmax=partsnap.size();
  HBTInt blocksize=nmax/nloop+1;
  for(HBTInt offset=0,iloop=0;offset<nmax;offset+=blocksize, iloop++)
  {
    HBTInt np=min(blocksize, nmax-offset);
    SnapshotOffset_t snapslice(offset, partsnap);//split the snapshot and do it piece by piece
    snapslice.ReSize(np);
    cout<<"i="<<iloop<<endl;
    
  #ifdef USE_LL
    SnapshotPos_t PartPos(snapslice);
    LinkedlistPara_t searcher(16, &PartPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
    cout<<"linked list compiled\n";
  #else
    GeoTree_t searcher;
    searcher.Build(snapslice);
    cout<<"tree built\n";
  #endif
    HBTInt iskip=partsnap.size()/SAMPLESIZE;
    #pragma omp parallel for
    for(int i=0;i<SAMPLESIZE;++i)
    {
      HBTInt pid=i*iskip;
      Profiles[i].Compute(partsnap.GetComovingPosition(pid), searcher);
    }
  }

  HBTInt NumProf[NBIN]={0};
  for(int pid=0;pid<SAMPLESIZE;++pid)
    for(int i=0;i<NBIN;i++)
      NumProf[i]+=Profiles[pid].n[i];
    
  float prof[NBIN];
  for(int i=0;i<NBIN;i++)
    prof[i]=NumProf[i]*1./SAMPLESIZE;
  
  save(prof, isnap);
    
  return 0;
}

void HaloSize_t::Compute(const HBTxyz &cen, LinkedlistPara_t &ll)
{
    vector <LocatedParticle_t> founds;
    founds.reserve(1024);
    ll.SearchSphereSerial(RMAX, cen, founds);
    for(auto &&p: founds)
    {
	int ibin=ceilf(logf(p.d2/RMIN/RMIN)/DlnX);
	if(ibin<0) ibin=0;
	else if(ibin>=NBIN) ibin=NBIN-1;
	n[ibin]++;
    }
}

void HaloSize_t::Compute(const HBTxyz &cen, GeoTree_t &tree)
{
    vector <LocatedParticle_t> founds;
    founds.reserve(1024);
    tree.Search(cen, RMAX, founds);
    for(auto &&p: founds)
    {
	int ibin=ceilf(logf(p.d2/RMIN/RMIN)/DlnX);
	if(ibin<0) ibin=0;
	else if(ibin>=NBIN) ibin=NBIN-1;
	n[ibin]++;
    }
}

void save(float prof[], int isnap)
{
  string filename=HBTConfig.SubhaloPath+"/postproc";
  mkdir(filename.c_str(), 0755);
  filename=filename+"/RandomCount_"+to_string(isnap)+".txt";
  ofstream file(filename);
  for(int i=0;i<NBIN;i++)
    file<<prof[i]<<" ";
  file<<endl;
  file.close();
}
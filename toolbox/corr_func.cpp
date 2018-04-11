/*to compute the average count profile around randomly selected particles
 */

#include <cmath>
#include <iostream>
#include <string>
#include <random>
#include <omp.h>

#include "../src/datatypes.h"
#include "../src/config_parser.h"
#include "../src/snapshot.h"
#include "../src/halo.h"
#include "../src/subhalo.h"
#include "../src/mymath.h"
#include "../src/linkedlist_parallel.h"
#include "../src/geometric_tree.h"

#define SAMPLESIZE 100000
// #define RMIN 1e-2 //only for deciding binwidth; not used to set the innermost bin edge.
// #define RMAX 2e1
// #define NBIN 20
#define RMIN 20e3 //only for deciding binwidth; not used to set the innermost bin edge.
#define RMAX 100e3
#define NBIN 10
#define USE_LL //algorithm: whether to use linkedlist or geotree for spatial search.
#define NLOOP 4
#define NDIV 128

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

class LogBinCollector_t: public ParticleCollector_t
{
  float DlnX;
  float Rmin2;
  int Nbin;
public:
  vector <HBTInt> N;
  LogBinCollector_t(float rmin, float rmax, int nbin)
  {
    Nbin=nbin;
    N.clear();
    N.resize(Nbin, 0);
    Rmin2=rmin*rmin;
    DlnX=logf(rmax/rmin)/(nbin-1)*2; //spacing in ln-space for r**2
  }
  void Collect(HBTInt index, HBTReal d2)
  {
    int ibin=ceilf(logf(d2/Rmin2)/DlnX);
    if(ibin<0) ibin=0;
    else if(ibin>=Nbin) ibin=Nbin-1;
    N[ibin]++;
  }
  void print(int ithread)
  {
    cout<<"ithread="<<ithread<<", Nbin="<<Nbin<<"N="<<N[0]<<","<<N[3]<<",rmin2="<<Rmin2<<endl;
  }
};

extern LogBinCollector_t Collector;
#pragma omp threadprivate(Collector)

LogBinCollector_t Collector(RMIN, RMAX, NBIN);

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
    
  omp_set_dynamic(0);//so that the threadprivate vars persist accross parallel regions
#pragma omp parallel
  Collector.print(omp_get_thread_num());
  
  const int nloop=NLOOP; //do it in blocks to save memory
  HBTInt nmax=partsnap.size();
  HBTInt blocksize=nmax/nloop+1;
  for(HBTInt offset=0,iloop=0;offset<nmax;offset+=blocksize, iloop++)
  {
    HBTInt np=min(blocksize, nmax-offset);
    SnapshotOffset_t snapslice(offset, partsnap);//split the snapshot and do it piece by piece
    snapslice.ReSize(np);
    cout<<"i="<<iloop<<endl;
    Timer_t timer;timer.Tick();
  #ifdef USE_LL
    SnapshotPos_t PartPos(snapslice);
    LinkedlistPara_t ll(NDIV, &PartPos, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);
    cout<<"linked list compiled\n";
  #else
    GeoTree_t tree;
    tree.Build(snapslice);
    cout<<"tree built\n";
  #endif
    timer.Tick();
    HBTInt iskip=partsnap.size()/SAMPLESIZE;
    #pragma omp parallel for schedule(dynamic, 10)
    for(int i=0;i<SAMPLESIZE;++i)
    {
      HBTInt pid=i*iskip;
      auto &cen=partsnap.GetComovingPosition(pid);
    #ifdef USE_LL
      ll.SearchSphereSerial(RMAX, cen, Collector);
    #else
      tree.Search(cen, RMAX, Collector);
    #endif
    }
    timer.Tick();cout<<"built in "<<timer.GetSeconds(1)<<" seconds, searched in "<<timer.GetSeconds(-1)<<" seconds\n";
  }

  HBTInt NumProf[NBIN]={0};
#pragma omp parallel
#pragma omp critical(reduce_count_profile)
    for(int i=0;i<NBIN;i++)
      NumProf[i]+=Collector.N[i];
    
  float prof[NBIN];
  for(int i=0;i<NBIN;i++)
    prof[i]=NumProf[i]*1./SAMPLESIZE;
  
  save(prof, isnap);
    
  return 0;
}

void save(float prof[], int isnap)
{
  string filename=HBTConfig.SubhaloPath+"/postproc";
  mkdir(filename.c_str(), 0755);
  filename=filename+"/RandomCount_"+to_string(isnap)+".outer.txt";
  ofstream file(filename);
  for(int i=0;i<NBIN;i++)
    file<<prof[i]<<" ";
  file<<endl;
  file.close();
}

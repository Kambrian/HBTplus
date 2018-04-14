#include <omp.h>
#include "linkedlist.h"

LinkedlistPara_t::LinkedlistPara_t(int ndiv, PositionData_t *data, HBTReal boxsize, bool periodic)
{
  #pragma omp parallel
  {
#ifdef _OPENMP
    int thread_id=omp_get_thread_num();
    int thread_num=omp_get_num_threads();
#else
    int thread_id=0, thread_num=1;
#endif
    #pragma omp single
    {
      LLs.resize(thread_num);
      Samples.resize(thread_num);
    }
    Samples[thread_id].init(thread_id, thread_num, data);
    LLs[thread_id].build(ndiv, &(Samples[thread_id]), boxsize, periodic);
  }
}
class SampleCollector_t: public ParticleCollector_t
{
  ParticleCollector_t &Collector;
  PositionSampleBase_t &Sample;
public:
  SampleCollector_t(ParticleCollector_t &collector, PositionSampleBase_t &sample):Collector(collector), Sample(sample)
  {}
  void Collect(HBTInt pid, HBTReal d2)
  {
    Collector.Collect(Sample.restore_id(pid), d2);
  }
};
void LinkedlistPara_t::SearchShellSerial(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter, ParticleCollector_t &collector)
{//serial version, which can be safely run inside another parallel region
  for(int thread_id=0;thread_id<LLs.size();thread_id++)
  {
    SampleCollector_t thread_collector(collector, Samples[thread_id]);
    LLs[thread_id].SearchShell(rmin, rmax, searchcenter, thread_collector);
  }
}

void Linkedlist_t::parallel_build(int ndiv, PositionData_t *data, HBTReal boxsize, bool periodic)
{
  assert(boxsize>0);//do not support auto-boxsize; has to make sure each thread use the same boxsize
  #pragma omp parallel
  {
#ifdef _OPENMP
    int thread_id=omp_get_thread_num();
    int thread_num=omp_get_num_threads();
#else
    int thread_id=0, thread_num=1;
#endif
    #pragma omp single
    {
      LLs.resize(thread_num);
      Samples.resize(thread_num);
    }
    Samples[thread_id].init(thread_id, thread_num, data);
    LLs[thread_id].build(ndiv, &(Samples[thread_id]), boxsize, periodic);
  }
  
  init(ndiv, data, boxsize, periodic);
  merge();
}
void Linkedlist_t::merge()
{//merge the LLs into the main list
  #pragma omp parallel for
  for(int ichain=0;ichain<HOC.size();ichain++)
  {
    auto &hoc=HOC[ichain];
    for(int ithread=0;ithread<LLs.size();ithread++)
    {
      auto pid=LLs[ithread].HOC[ichain];
      while(pid>=0)
      {
	{//push back
	auto true_pid=Samples[ithread].restore_id(pid);
	List[true_pid]=hoc;
	hoc=true_pid;
	}
	pid=LLs[ithread].List[pid];//next particle
      }
    }
  }
  
  //equivalent way:
//   #pragma omp parallel for
//   for(int ichain=0;ichain<HOC.size();ichain++)
//   {
//     auto *p=&HOC[ichain];
//     for(int ithread=0;ithread<LLs.size();ithread++)
//     {
//       auto pid=LLs[ithread].HOC[ichain];
//       while(pid>=0)//copy till tail
//       {
// 	*p=Samples[ithread].restore_id(pid);//copy
// 	pid=LLs[ithread].List[pid];//next value
// 	p=&List[*p];//next storage
//       }
//     }
//     *p=-1;//close the chain
//   }

  LLs.clear();
  Samples.clear();
}

void LinkedlistLinkGroup(HBTReal radius, const Snapshot_t &snapshot, vector <HBTInt> &GrpLen, vector <HBTInt> &GrpTags, int ndiv)
/* link particles in the given snapshot into groups.
 * Output: filled GrpLen and GrpTags (0~Ngroups-1), down to mass=1 (diffuse particles)
 * */
{
GrpTags.assign(snapshot.size(), -1);

cout<<"Building linkedlist...\n"<<flush;

SnapshotPos_t posdata(snapshot);
Linkedlist_t ll(ndiv, &posdata, HBTConfig.BoxSize, HBTConfig.PeriodicBoundaryOn);

cout<<"Linking Groups...\n"<<flush;
HBTInt grpid=0;
HBTInt printstep=snapshot.size()/100, progress=printstep;
cout<<"00%"<<flush;
for(HBTInt i=0;i<snapshot.size();i++)
{
	if(GrpTags[i]<0)
	{
		if(i>=progress)
		{
		  cout<<"\b\b\b"<<setw(2)<<progress/printstep<<"%"<<flush;
		  progress+=printstep;
		}
		GrpTags[i]=grpid; //infect the seed particle
		HBTInt grplen=1+ll.TagFriendsOfFriends(i,grpid, GrpTags, radius);
		GrpLen.push_back(grplen);
		grpid++;
	}
}

cout<<"Found "<<GrpLen.size()<<" Groups\n";
}

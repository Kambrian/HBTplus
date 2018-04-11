#include "mymath.h"
#include "linkedlist.h"

class PositionSample_t: public PositionData_t
{
private:
  int ThreadId, NumThreads;
  PositionData_t *Data;
  HBTInt np;
public:
  void init(int ithread, int nthread,  PositionData_t *data)
  {
    Data=data;
    ThreadId=ithread;
    NumThreads=nthread;
    HBTInt n0=Data->size();
    np=n0/nthread+((n0%nthread)>ithread);
  }
  const HBTxyz & operator [](HBTInt i) const
  {
    return (*Data)[restore_id(i)];
  }
  size_t size() const
  {
    return np;
  }
  HBTInt restore_id(HBTInt i) const
  {
    return i*NumThreads+ThreadId;
  }
  void restore_id(vector <LocatedParticle_t> &particles) const
  {
    for(auto &&p: particles)
      p.index=restore_id(p.index);
  }
};

class LinkedlistPara_t
{
private:
  vector<Linkedlist_t> LLs;
  vector <PositionSample_t> Samples;
public:
  LinkedlistPara_t(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false); 
  template <class REDUCIBLECOLLECTOR>
  void SearchSphere(HBTReal radius, const HBTxyz &searchcenter, ParticleCollector_t &collector)
  {
    SearchShell<REDUCIBLECOLLECTOR>(-1., radius, searchcenter, collector);
  }
  void SearchSphereSerial(HBTReal radius, const HBTxyz &searchcenter, ParticleCollector_t &collector)
  {
    SearchShellSerial(-1., radius, searchcenter, collector);
  }
  template <class REDUCIBLECOLLECTOR>
  void SearchShell(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter, ParticleCollector_t &collector);
  void SearchShellSerial(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter, ParticleCollector_t &collector);
};

class ReducibleCollector_t: public ParticleCollector_t
{
  virtual void Reduce(ParticleCollector_t &final_collector)=0;//defines how to add the results from each thread together into final_collector.
};

template <class REDUCIBLECOLLECTOR>
class ReducibleSampleCollector_t: public ParticleCollector_t
{
  REDUCIBLECOLLECTOR Collector;
  PositionSample_t &Sample;
public:
  ReducibleSampleCollector_t(PositionSample_t &sample): Sample(sample), Collector()
  {}
  void Collect(HBTInt id, HBTReal d2)
  {
    Collector.Collect(Sample.restore_id(id), d2);
  }
  void Reduce(ParticleCollector_t &collector)
  {
    Collector.Reduce(collector);
  }
};
  
template <class REDUCIBLECOLLECTOR>
void LinkedlistPara_t::SearchShell(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter, ParticleCollector_t &collector)
{//parallel version. not suitable for use inside another parallel region. must specify a ReducibleCollector_t template parameter to define Collect and Reduce method for each thread collector.
#pragma omp parallel for
  for(int thread_id=0;thread_id<LLs.size();thread_id++)
  {
    ReducibleSampleCollector_t<REDUCIBLECOLLECTOR> thread_founds(Samples[thread_id]);
    LLs[thread_id].SearchShell(rmin, rmax, searchcenter, thread_founds);
    #pragma omp critical(insert_linklist_founds) //this prevents nested parallelization
    {
      thread_founds.Reduce(collector);
    }
  }
}

class ReducibleLocatedParticleCollector_t:  public LocatedParticleCollector_t, public ReducibleCollector_t
/* a simple collector to be used for parallel search; it keeps a local vector to store located particles in each thread, and then dump them to the final output collector*/
{
public:
  ReducibleLocatedParticleCollector_t(HBTInt n_reserve = 0):LocatedParticleCollector_t(n_reserve)
  {}
  void Collect(HBTInt index, HBTReal d2)
  {
    LocatedParticleCollector_t::Collect(index, d2);
  }
  void Reduce(ParticleCollector_t &collector)
  {
    for(auto &&p: Founds)
      collector.Collect(p.index, p.d2);
  }
};



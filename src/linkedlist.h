#ifndef LINKEDLIST_HEADER_INCLUDED
#define LINKEDLIST_HEADER_INCLUDED
#include "mymath.h"
#include "linkedlist_base.h"
class PositionSampleBase_t: public PositionData_t
{// a sample created from PositionData_t
protected:
  PositionData_t *Data;
  HBTInt np;
public:
  const HBTxyz & operator [](HBTInt i) const
  {
    return (*Data)[restore_id(i)];
  }
  size_t size() const
  {
    return np;
  }
  virtual HBTInt restore_id(HBTInt i) const//has to make it virtual to be overriden by derived class. otherwise the derived class will use the wrong [].
  {
    return i;
  }
  void restore_id(vector <LocatedParticle_t> &particles) const
  {
    for(auto &&p: particles)
      p.index=restore_id(p.index);
  }
};
class PositionSampleLattice_t: public PositionSampleBase_t
{//sample by skipping. suitable for searching multiple samples in parallel
  int ThreadId, NumThreads;
public:
  void init(int ithread, int nthread,  PositionData_t *data)
  {
    Data=data;
    ThreadId=ithread;
    NumThreads=nthread;
    HBTInt n0=Data->size();
    np=n0/nthread+((n0%nthread)>ithread);
  }
  HBTInt restore_id(HBTInt i) const//has to make it virtual to be overriden by derived class. otherwise the derived class will use the wrong [].
  {
    return i*NumThreads+ThreadId;
  }
};
class PositionSampleBlock_t: public PositionSampleBase_t
{//sample in blocks. suitable for merging.
  HBTInt offset;
public:
  void init(int ithread, int nthread,  PositionData_t *data)
  {
    AssignTasks(ithread, nthread, data->size(), offset, np);
    np-=offset;
    Data=data;
  }
  HBTInt restore_id(HBTInt i) const
  {
    return offset+i;
  }
};

class LinkedlistPara_t
{//built and searched in parallel. can be searched in serial as well. lower efficiency than Linkedlist_t, especially when built with larger number of threads.
private:
  vector<LinkedlistBase_t> LLs;
  vector <PositionSampleLattice_t> Samples;
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
  template <class REDUCIBLECOLLECTOR>
  void SearchCylinder(HBTReal radius_z, HBTReal radius_p, const HBTxyz &searchcenter, ParticleCollector_t &collector);
  void SearchCylinderSerial(HBTReal radius_z, HBTReal radius_p, const HBTxyz &searchcenter, ParticleCollector_t &collector);
};

class ReducibleCollector_t: public ParticleCollector_t
{
  virtual void Reduce(ParticleCollector_t &final_collector)=0;//defines how to add the results from each thread together into final_collector.
};

template <class REDUCIBLECOLLECTOR>
class ReducibleSampleCollector_t: public ParticleCollector_t
{
  REDUCIBLECOLLECTOR Collector;
  PositionSampleBase_t &Sample;
public:
  ReducibleSampleCollector_t(PositionSampleBase_t &sample): Sample(sample), Collector()
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

template <class REDUCIBLECOLLECTOR>
void LinkedlistPara_t::SearchCylinder(HBTReal radius_z, HBTReal radius_p, const HBTxyz &searchcenter, ParticleCollector_t &collector)
{//parallel version. not suitable for use inside another parallel region. must specify a ReducibleCollector_t template parameter to define Collect and Reduce method for each thread collector.
#pragma omp parallel for
  for(int thread_id=0;thread_id<LLs.size();thread_id++)
  {
    ReducibleSampleCollector_t<REDUCIBLECOLLECTOR> thread_founds(Samples[thread_id]);
    LLs[thread_id].SearchCylinder(radius_z, radius_p, searchcenter, thread_founds);
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

class Linkedlist_t:public LinkedlistBase_t
{//built in parallel
private:
  vector<LinkedlistBase_t> LLs;
  vector <PositionSampleBlock_t> Samples;
  void merge();
public:
  Linkedlist_t():LinkedlistBase_t()
  {}
  Linkedlist_t(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false, bool build_in_parallel=true):LinkedlistBase_t()
  {
    if(build_in_parallel)
      parallel_build(ndiv, data, boxsize, periodic);
    else
      build(ndiv, data, boxsize, periodic);
  }
  void parallel_build(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false);
};

extern void LinkedlistLinkGroup(HBTReal radius, const Snapshot_t &snapshot, vector <HBTInt> &GrpLen, vector <HBTInt> &GrpTags, int ndiv=256);
#endif

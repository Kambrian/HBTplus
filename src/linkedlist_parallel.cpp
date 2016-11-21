#include "mymath.h"
#include "linkedlist.cpp"

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
    HBTInt n0=Data.size();
    np=n0/nthread+((n0%nthread)>ithread);
  }
  const HBTxyz & operator [](HBTInt i) const
  {
    return (*Data)[i*NumThreads+ThreadId];
  }
  size_t size() const
  {
    return np;
  }
};

Linkedlist_t ThreadLL;
PositionSample_t ThreadSample;
#pragma omp thread_private(ThreadLL, ThreadSample)

class LinkedlistPara_t
{
private:
  LinkedlistPara_t(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false): 
  {
    #pragma omp parallel
    {
      int thread_id=omp_get_thread_num();
      int thread_num=omp_get_num_threads();
      ThreadSample.init(thread_id, thread_num, data);
      ThreadLL.init(ndiv, &ThreadSample, boxsize, periodic);
    }
  }
  void SearchSphere(HBTReal radius, const HBTxyz &searchcenter, vector <HBTInt> &found_ids, int nmax_guess=8)
  {
    found_ids.clear();
#pragma omp parallel
    {
      vector <HBTInt> thread_founds(nmax_guess);
      ThreadLL.SearchSphere(radius, searchcenter, thread_founds);
      #pragma omp critical
      {
	found_ids.insert(found_ids.end(), thread_founds.begin(), thread_founds.end());
      }
    }
  }
};



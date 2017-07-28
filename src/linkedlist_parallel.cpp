#include <omp.h>
#include "linkedlist_parallel.h"

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
void LinkedlistPara_t::SearchShell(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter, vector <LocatedParticle_t> &founds)
{//parallel version. not suitable for use inside another parallel region.
#pragma omp parallel for
  for(int thread_id=0;thread_id<LLs.size();thread_id++)
  {
    vector <LocatedParticle_t> thread_founds;
    LLs[thread_id].SearchShell(rmin, rmax, searchcenter, thread_founds);
    Samples[thread_id].restore_id(thread_founds);
    #pragma omp critical(insert_linklist_founds) //this prevents nested parallelization
    {
      founds.insert(founds.end(), thread_founds.begin(), thread_founds.end());
    }
  }
}
void LinkedlistPara_t::SearchShellSerial(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter, vector <LocatedParticle_t> &founds)
{//serial version, which can be safely run inside another parallel region
  for(int thread_id=0;thread_id<LLs.size();thread_id++)
  {
    vector <LocatedParticle_t> thread_founds;
    LLs[thread_id].SearchShell(rmin, rmax, searchcenter, thread_founds);
    Samples[thread_id].restore_id(thread_founds);
    founds.insert(founds.end(), thread_founds.begin(), thread_founds.end());
  }
}


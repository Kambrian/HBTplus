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
class SampleCollector_t: public ParticleCollector_t
{
  ParticleCollector_t &Collector;
  PositionSample_t &Sample;
public:
  SampleCollector_t(ParticleCollector_t &collector, PositionSample_t &sample):Collector(collector), Sample(sample)
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
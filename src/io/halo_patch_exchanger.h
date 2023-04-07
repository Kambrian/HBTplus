/* common functions for exchanging and combining halo segments read in parallel on different processors
*/
#ifndef HALO_PATCH_EXCHANGER_INCLUDED
#define HALO_PATCH_EXCHANGER_INCLUDED
#include "../hdf_wrapper.h"
#include "../halo.h"
#include "../mpi_wrapper.h"
#include "../halo_particle_iterator.h"

namespace HaloPatchExchanger
{
  //this func differs from the particle_exchanger one in that the halos are fragmented onto different processors here.
  extern void ExchangeAndMerge(MpiWorker_t& world, vector< Halo_t >& Halos);
}
#endif

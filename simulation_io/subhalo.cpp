#include <iostream>
#include <new>

#include "../datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"

void SubHaloSnapshot_t::descend_into(HaloSnapshot_t &halo_snap, Snapshot_t &part_snap)
{
  vector<HBTInt> ParticleToHost(part_snap.GetNumberOfParticles(), SpecialConst::NullHaloId);
  for(int haloid=0;haloid<halo_snap.Halos.Size();haloid++)
  {
	for(int i=0;i<halo_snap.Halos[haloid].Particles.Size();i++)
	  ParticleToHost[i]=haloid;
  }
  for(int subid=0;subid<SubHalos.Size();subid++)
  {
	//rely on track_particle
	SubHalos[subid].HostHaloId=ParticleToHost[SubHalos[subid].TrackParticleId];
	//alternatives: CoreTrack; Split;
  }
  //alternative: trim particles outside fof
}

void SubHaloSnapshot_t::decide_centrals()
{
  for(HBTInt subid=0;subid<SubHalos.Size();subid++)
  {//FIXME.............
  }
}

void SubHaloSnapshot_t::refine_particles()
{//it's more expensive to build an exclusive list. so do inclusive here. 
  //TODO: ensure the inclusive unbinding is stable (contaminating particles from big subhaloes may hurdle the unbinding
  for(int subid=0;subid<SubHalos.Size();subid++)
  {
	SubHalos[subid].unbind();
  }
}

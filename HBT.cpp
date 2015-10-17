using namespace std;
#include <iostream>
#include <string>
#include <cstdlib>

#include "simulation_io/snapshot.h"
#include "simulation_io/halo.h"
#include "simulation_io/subhalo.h"
#include "mymath.h"

int main(int argc, char **argv)
{
#ifdef _OPENMP
 omp_set_nested(0);
#endif
  int snapshot_start, snapshot_end;
  ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  mkdir(HBTConfig.SubhaloPath.c_str(), 0755);
  MarkHBTVersion();
	
  SubHaloSnapshot_t subsnap;
  
  subsnap.Load(snapshot_start-1, true);
    
  for(int isnap=snapshot_start;isnap<=snapshot_end;isnap++)
  {
	Snapshot_t partsnap;
	partsnap.Load(isnap);
	subsnap.SetSnapshotIndex(isnap);
	
	HaloSnapshot_t halosnap;
	halosnap.Load(isnap);
	#pragma omp parallel
	{
	halosnap.ParticleIdToIndex(partsnap);
	subsnap.ParticleIdToIndex(partsnap);
	subsnap.AssignHosts(halosnap);
	subsnap.PrepareCentrals(halosnap);
	}
	
	subsnap.RefineParticles();
	
	#pragma omp parallel
	{
	subsnap.UpdateTracks();
	subsnap.ParticleIndexToId();
	}
	
	subsnap.Save();
  }
  
  return 0;
}
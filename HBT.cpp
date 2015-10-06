using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "simulation_io/snapshot.h"
#include "simulation_io/halo.h"
#include "simulation_io/subhalo.h"
#include "mymath.h"

int main(int argc, char **argv)
{
  int snapshot_start, snapshot_end;
  ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  SubHaloSnapshot_t subsnap;
  
  subsnap.Load(snapshot_start-1);
    
  for(int isnap=snapshot_start;isnap<=snapshot_end;isnap++)
  {
	Snapshot_t partsnap;
	HaloSnapshot_t halosnap;
	partsnap.Load(isnap);
	halosnap.Load(isnap);
	halosnap.ParticleIdToIndex(partsnap);
	subsnap.ParticleIdToIndex(partsnap);
	
	subsnap.AssignHost(halosnap);
	
	halosnap.AverageCoordinates();
	subsnap.AverageCoordinates();
	subsnap.DecideCentrals(halosnap);
	subsnap.FeedCentrals(halosnap);
	subsnap.RefineParticles();
	
	subsnap.ExtendTrackIds();
	
	subsnap.ParticleIndexToId();
	subsnap.Save();
  }
  
  return 0;
}

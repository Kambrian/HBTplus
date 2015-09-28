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
  
  subsnap.Load(HBTConfig, snapshot_start-1);
  
  for(int isnap=snapshot_start;isnap<=snapshot_end;isnap++)
  {
	Snapshot_t partsnap;
	HaloSnapshot_t halosnap;
	partsnap.Load(HBTConfig, isnap);
	halosnap.Load(HBTConfig, isnap);
	halosnap.ParticleIdToIndex(partsnap);
	subsnap.ParticleIdToIndex(partsnap);
		
	subsnap.descend_into(halosnap, partsnap);
	
	halosnap.AverageHaloCoordinates(partsnap);
	subsnap.decide_centrals();
	
	subsnap.refine_particles();
	
	subsnap.ParticleIndexToId(partsnap);
	subsnap.Save(HBTConfig);
// 	subsnap_old=subsnap;
  }
  
  return 0;
}

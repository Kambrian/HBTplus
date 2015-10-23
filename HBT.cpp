using namespace std;
#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>

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
	
  SubhaloSnapshot_t subsnap;
  
  subsnap.Load(snapshot_start-1, true);
  
  vector <chrono::high_resolution_clock::time_point> tickers(20);
  int ticker_count=0;
  #define TICK_TIME() tickers[ticker_count++]=chrono::high_resolution_clock::now()
  ofstream time_log(HBTConfig.SubhaloPath+"/timing.log", fstream::out|fstream::app);
  
  for(int isnap=snapshot_start;isnap<=snapshot_end;isnap++)
  {
	TICK_TIME();
	ParticleSnapshot_t partsnap;
	partsnap.Load(isnap);
	subsnap.SetSnapshotIndex(isnap);
	
	HaloSnapshot_t halosnap;
	halosnap.Load(isnap);
	TICK_TIME();
	#pragma omp parallel
	{
	halosnap.ParticleIdToIndex(partsnap);
	subsnap.ParticleIdToIndex(partsnap);
	#pragma omp master
	TICK_TIME();
	subsnap.AssignHosts(halosnap);
	subsnap.PrepareCentrals(halosnap);
	}

	TICK_TIME();
	subsnap.RefineParticles();
	TICK_TIME();
	
	#pragma omp parallel
	{
	subsnap.UpdateTracks();
	subsnap.ParticleIndexToId();
	}
	TICK_TIME();
	
	subsnap.Save();
	
	TICK_TIME();
	time_log<<isnap;
	for(int i=1;i<ticker_count;i++)
	  time_log<<"\t"<<chrono::duration_cast<chrono::duration<double> >(tickers[i]-tickers[i-1]).count();
	time_log<<endl;
	ticker_count=0;
  }
  
  return 0;
}
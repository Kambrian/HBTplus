using namespace std;
#include <cstdlib>
#include <iostream>
#include <string>
#include <omp.h>

#include "src/mpi_wrapper.h"
#include "src/datatypes.h"
#include "src/config_parser.h"
#include "src/snapshot.h"
#include "src/halo.h"
#include "src/subhalo.h"
#include "src/mymath.h"
#include "src/particle_exchanger.h"

int main(int argc, char **argv)
{
 MPI_Init(&argc, &argv);
 MpiWorker_t world(MPI_COMM_WORLD);
#ifdef _OPENMP
 omp_set_nested(0);
#endif
   
  int snapshot_start, snapshot_end;
  if(0==world.rank())
  {
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
	mkdir(HBTConfig.SubhaloPath.c_str(), 0755);
	HBTConfig.DumpParameters();
  }
  HBTConfig.BroadCast(world, 0, snapshot_start, snapshot_end);

  SubhaloSnapshot_t subsnap;
  
  subsnap.Load(world, snapshot_start-1, true);
  
  Timer_t timer;
  ofstream time_log;
  if(world.rank()==0)
  {
  time_log.open(HBTConfig.SubhaloPath+"/timing.log", fstream::out|fstream::app);
  time_log<<fixed<<setprecision(1);//<<setw(8);
  }
  for(int isnap=snapshot_start;isnap<=snapshot_end;isnap++)
  {
	timer.Tick(world.Communicator);
	ParticleSnapshot_t partsnap;
	partsnap.Load(world, isnap);
	subsnap.SetSnapshotIndex(isnap);
	HaloSnapshot_t halosnap;
	halosnap.Load(world, isnap);
	
	timer.Tick(world.Communicator);
// 	cout<<"updating halo particles...\n";
	halosnap.UpdateParticles(world, partsnap);
	timer.Tick(world.Communicator);
// 	if(world.rank()==0) cout<<"updateing subsnap particles...\n";
	subsnap.UpdateParticles(world, partsnap);
	
	timer.Tick(world.Communicator);
	subsnap.AssignHosts(world, halosnap, partsnap);
	subsnap.PrepareCentrals(halosnap);
/*
	timer.Tick(world.Communicator);
	subsnap.RefineParticles();
	
	timer.Tick(world.Communicator);
	subsnap.UpdateTracks(world, halosnap);
	
	cout<<" start to save "<<subsnap.Subhalos.size()<<" subs on thread "<<world.rank()<<endl;
	
	timer.Tick(world.Communicator);
	subsnap.Save(world);
	*/
	timer.Tick(world.Communicator);
	if(world.rank()==0)
  {
	time_log<<isnap;
	for(int i=1;i<timer.Size();i++)
	  time_log<<"\t"<<timer.GetSeconds(i);
	time_log<<endl;
	timer.Reset();
  }
  }
  
  MPI_Finalize();
  return 0;
}
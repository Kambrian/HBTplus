using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "../src/snapshot.h"
#include "../src/halo.h"
#include "../src/mymath.h"

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

    cout<<argv[0]<<" run using "<<world.size()<<" mpi tasks";
    #ifdef _OPENMP
    #pragma omp parallel
    #pragma omp master
    cout<<", each with "<<omp_get_num_threads()<<" threads";
    #endif
    cout<<endl;
  }
  HBTConfig.BroadCast(world, 0, snapshot_start, snapshot_end);

  int isnap=snapshot_start;

  Timer_t timer;
  ParticleSnapshot_t partsnap;
  HaloSnapshot_t halosnap;

  timer.Tick(world.Communicator);
  partsnap.Load(world, isnap, false);
  timer.Tick(world.Communicator);
  if(world.rank()==0)
    cout<<"\nTiming:>>>\n\tSnapshot loaded in "<<timer.GetSeconds()<<" secs\n\n";
  stringstream ss;
  ss<<"rank "<<world.rank()<<": "<<partsnap.Particles.size()<<" particles, 0,100,end are: "<<partsnap.Particles[0]<<endl<<partsnap.Particles[100]<<endl<<partsnap.Particles.back()<<endl;
  cout<<ss.str()<<flush;
  MPI_Barrier(world.Communicator);

  timer.Tick(world.Communicator);
  halosnap.Load(world, partsnap);
  timer.Tick(world.Communicator);
  if(world.rank()==0)
    cout<<"\nTiming:>>>\n\tHalos loaded in "<<timer.GetSeconds()<<" secs\n\n";

  ss<<"rank "<<world.rank()<<": "<<halosnap.Halos.size()<<" halos. \n";
  if(halosnap.Halos.size()>0)
  ss<<"Halo 0: Id="<<halosnap.Halos[0].HaloId<<", "<<halosnap.Halos[0].Particles.size()<<" particles: "<<halosnap.Halos[0].Particles[0]<<"; "<<halosnap.Halos[0].Particles[1]<<";... "<<halosnap.Halos[0].Particles.back()<<endl;
  if(halosnap.Halos.size()>1)
  ss<<"Halo 1: Id="<<halosnap.Halos[1].HaloId<<", "<<halosnap.Halos[1].Particles.size()<<" particles: "<<halosnap.Halos[1].Particles[0]<<"; "<<halosnap.Halos[1].Particles[1]<<";... "<<halosnap.Halos[1].Particles.back()<<endl;
  if(halosnap.Halos.size()>1)
  ss<<"Halo -1: Id="<<halosnap.Halos.back().HaloId<<", "<<halosnap.Halos.back().Particles.size()<<" particles: "<<halosnap.Halos.back().Particles[0]<<"; "<<halosnap.Halos.back().Particles[1]<<";... "<<halosnap.Halos[1].Particles.back()<<endl;
  cout<<ss.str()<<flush;
  MPI_Barrier(world.Communicator);

  halosnap.Clear();
  HBTConfig.GroupFileFormat="gadget4_hdf2";
  timer.Tick(world.Communicator);
  halosnap.Load(world, partsnap);
  timer.Tick(world.Communicator);
  if(world.rank()==0)
    cout<<"\nTiming:>>>\n\t halosnap loaded with hdf2 in "<<timer.GetSeconds()<<" secs\n\n";

  ss<<"rank "<<world.rank()<<": "<<halosnap.Halos.size()<<" halos. \n";
  if(halosnap.Halos.size()>0)
  ss<<"Halo 0: Id="<<halosnap.Halos[0].HaloId<<", "<<halosnap.Halos[0].Particles.size()<<" particles: "<<halosnap.Halos[0].Particles[0]<<"; "<<halosnap.Halos[0].Particles[1]<<";... "<<halosnap.Halos[0].Particles.back()<<endl;
  if(halosnap.Halos.size()>1)
  ss<<"Halo 1: Id="<<halosnap.Halos[1].HaloId<<", "<<halosnap.Halos[1].Particles.size()<<" particles: "<<halosnap.Halos[1].Particles[0]<<"; "<<halosnap.Halos[1].Particles[1]<<";... "<<halosnap.Halos[1].Particles.back()<<endl;
  if(halosnap.Halos.size()>1)
  ss<<"Halo -1: Id="<<halosnap.Halos.back().HaloId<<", "<<halosnap.Halos.back().Particles.size()<<" particles: "<<halosnap.Halos.back().Particles[0]<<"; "<<halosnap.Halos.back().Particles[1]<<";... "<<halosnap.Halos[1].Particles.back()<<endl;
  cout<<ss.str();

  return 0;
}

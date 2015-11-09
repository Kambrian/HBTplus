using namespace std;
#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;
// #include "mpi.h"

#include "src/snapshot.h"
#include "src/halo.h"
#include "src/subhalo.h"
#include "src/mymath.h"

int main(int argc, char **argv)
{
  mpi::environment env;
  mpi::communicator world;
#ifdef _OPENMP
 omp_set_nested(0);
#endif/*
  int  numtasks, rank, len, rc; 
  char hostname[MPI_MAX_PROCESSOR_NAME];
  rc = MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Get_processor_name(hostname, &len);*/
   
  int snapshot_start, snapshot_end;
  if(0==world.rank())
  {
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
	mkdir(HBTConfig.SubhaloPath.c_str(), 0755);
	MarkHBTVersion();
  }
  broadcast(world, HBTConfig, 0);
  
  cout<< HBTConfig.SnapshotPath<< " from "<<world.rank()<<" of "<<world.size()<<" on "<<env.processor_name()<<endl;
  cout<< HBTConfig.SnapshotIdList[10]<< " from "<<world.rank()<<" of "<<world.size()<<" on "<<env.processor_name()<<endl;
  cout<< HBTConfig.IsSet[2]<< " from "<<world.rank()<<" of "<<world.size()<<" on "<<env.processor_name()<<endl;
  
  return 0;
}
using namespace std;
#include <iostream>
#include <string>
#include <cstdlib>
#include <omp.h>
#include "mpi.h"

#include "src/datatypes.h"
#include "src/config_parser.h"
#include "src/boost_mpi.h"
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
#endif
   /*
  int snapshot_start, snapshot_end;
  if(0==world.rank())
  {
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
	mkdir(HBTConfig.SubhaloPath.c_str(), 0755);
	MarkHBTVersion();
  }
  broadcast(world, HBTConfig, 0);
  
  cout<< HBTConfig.SnapshotPath<< " from "<<world.rank()<<" of "<<world.size()<<" on "<<env.processor_name()<<endl;
  cout<< HBTConfig.SnapshotIdList<< " from "<<world.rank()<<" of "<<world.size()<<" on "<<env.processor_name()<<endl;
  cout<< HBTConfig.IsSet[2]<< " from "<<world.rank()<<" of "<<world.size()<<" on "<<env.processor_name()<<endl;
  
*/
  
  typedef vector <int> ParticleList_t;
  
#define MSG_LEN 100000
  ParticleList_t sendbuf(MSG_LEN, 1), recvbuf(MSG_LEN);
  /*
  if(world.rank()==0)
	MPI_Send(sendbuf.data(), MSG_LEN, MPI_INT, 1, 0, MPI_COMM_WORLD);
  else if(world.rank()==1)
	MPI_Recv(recvbuf.data(), MSG_LEN, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
  */
  if(world.rank()==0)
 	world.send(1, 0, sendbuf);
  else if(world.rank()==1)
	world.recv(0, 0, recvbuf);
  
  cout<<recvbuf[0]<<" on "<<world.rank()<<endl;
  /*
  vector <ParticleList_t> SendCells(world.size(), ParticleList_t());
  for(int i=0;i<SendCells.size();i++)
  SendCells[i].assign(MSG_LEN/world.size(), world.rank());
  cout<<"Ready to send on "<<world.rank()<<endl;
  vector <int> SendSizes, ReceiveSizes(world.size());
  for(auto && c :SendCells)
	SendSizes.push_back(c.size());
//  MPI_Alltoall(SendSizes.data(),
  
  vector <ParticleList_t> ReceiveCells(world.size(), ParticleList_t());
  for(int i=0;i<ReceiveCells.size();i++)
	ReceiveCells.resize(ReceiveSizes[i]);
  
//   all_to_all(world, SendCells, ReceiveCells);
//   SendCells.clear();
  
  for(int i=0;i<world.size();i++)
	cout<<world.rank()<<" : "<<ReceiveCells[i].size()<<endl;
    */
  return 0;
}
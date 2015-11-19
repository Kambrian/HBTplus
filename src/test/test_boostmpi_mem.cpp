using namespace std;
#include <iostream>
#include <string>
#include "mpi.h"

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;

int main(int argc, char **argv)
{
  mpi::environment env;
  mpi::communicator world;
    
#define MSG_LEN 100000
  vector <int> sendbuf(MSG_LEN, 1), recvbuf(MSG_LEN);

#define USE_BOOST
  
#ifndef USE_BOOST
  if(world.rank()==0)
	MPI_Send(sendbuf.data(), MSG_LEN, MPI_INT, 1, 0, MPI_COMM_WORLD);
  else if(world.rank()==1)
	MPI_Recv(recvbuf.data(), MSG_LEN, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
#else
  if(world.rank()==0)
	world.send(1, 0, sendbuf);
  else if(world.rank()==1)
	world.recv(0, 0, recvbuf);
#endif
  
  if(world.rank()==1)
	cout<<"Data received: "<<recvbuf[0]<<","<<recvbuf[1]<<"...\n";

  return 0;
}
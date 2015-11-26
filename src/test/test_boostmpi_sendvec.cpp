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
    
#define MSG_LEN 100
  vector <int> sendbuf(MSG_LEN, 1), recvbuf(MSG_LEN);

  if(world.rank()==0)
	world.send(1, 0, sendbuf.data(), sendbuf.size());
  else if(world.rank()==1)
	//you cannot send array and receive vec. have to match types: world.recv(0, 0, recvbuf) will fail.
	world.recv(0, 0, recvbuf.data(), recvbuf.size());
  
  if(world.rank()==1)
	cout<<recvbuf.size()<<" elements received: "<<recvbuf[0]<<","<<recvbuf[1]<<"...\n";

  return 0;
}
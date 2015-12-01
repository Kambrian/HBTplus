using namespace std;
#include <iostream>
#include <string>
// #include "mpi.h"

#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

int main(int argc, char **argv)
{
  mpi::environment env;
  mpi::communicator world;
    
#define MSG_LEN 1000000
  vector <int> sendbuf(MSG_LEN, 1), recvbuf(MSG_LEN);

  MPI_Comm comm=world;
  if(world.rank()==0)
  {
	mpi::packed_oarchive oa(comm);
	oa << sendbuf;
	auto sendptr = const_cast<void*>(oa.address());
	// cast to int because MPI uses ints for sizes like it's still 1990
	int sendsize = static_cast<int>(oa.size());
	MPI_Send(&sendsize, 1, MPI_INT, 1, 0, comm);
	MPI_Ssend(sendptr, sendsize, MPI_PACKED, 1, 0, comm);
	cout<<sendsize<<","<<oa.size()<<endl;
  }
  else if (world.rank()==1)
  {
	mpi::packed_iarchive ia(comm);
	int recvsize;
	MPI_Recv(&recvsize, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
	ia.resize(recvsize);
	auto recvptr = ia.address();
	MPI_Recv(recvptr, recvsize, MPI_PACKED, 0, 0, comm, MPI_STATUS_IGNORE);
	ia >> recvbuf;
	cout<<"Data received: "<<recvbuf[0]<<","<<recvbuf[1]<<"...\n";
  }
  
  return 0;
}
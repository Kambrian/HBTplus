using namespace std;
#include <iostream>
#include <string>
#include "mpi.h"
#include <iterator>
#include <vector>

int main(int argc, char **argv)
{
  int myrank, x;
  MPI_Init(&argc, &argv);  
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
  
  if(myrank == 0) 
  { 
	x=1000;
  } 
  MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
  cout<<"x= "<<x<<" on "<<myrank<<endl;

  MPI_Request Req;
  const int nmax=100;
  if(myrank==0)
  {
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0;i<nmax;i++)
	MPI_Send(&i, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	
  }
  else
  {
	vector <int> y(nmax);
	vector <MPI_Request> Reqs(nmax);
	for(int i=0;i<nmax;i++)
	{
	  MPI_Irecv(&y[i], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Reqs[i]);
	  cout<<y[i]<<" ";
	}
	cout<<endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Waitall(Reqs.size(), Reqs.data(), MPI_STATUSES_IGNORE);//without this, the receive may not happen inside this block, where y is destroyed.
	copy(y.cbegin(), y.cend(), ostream_iterator<int>(cout, ", "));
	cout<<endl;
  }
  
	MPI_Finalize();
	return 0;
}
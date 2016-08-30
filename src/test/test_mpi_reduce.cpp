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
  
  x=myrank;
  if(myrank==0)
    MPI_Reduce(MPI_IN_PLACE, &x, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce(&x, &x, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
  cout<<"thread "<<myrank<<": "<<x<<endl;
  
  MPI_Finalize();
  return 0;
}
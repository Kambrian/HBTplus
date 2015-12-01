using namespace std;
#include <iostream>
#include <string>
#include "mpi.h"

int main(int argc, char **argv)
{
#define MSG_LEN 1000000
#define BUF_LEN (4*(MSG_LEN+2))
  int position, i, j[MSG_LEN], a[2]; 
  char buff[BUF_LEN]; 
  int myrank;
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
  if(myrank == 0) 
  { 
	position = 0; 
	MPI_Pack(&a, 2, MPI_INT, buff, BUF_LEN, &position, MPI_COMM_WORLD); 
// 	cout<<position<<endl;
	MPI_Pack(j, MSG_LEN, MPI_INT, buff, BUF_LEN, &position, MPI_COMM_WORLD); 
// 	cout<<position<<endl;
	MPI_Send(&position, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	MPI_Send( buff, position, MPI_PACKED, 1, 0, MPI_COMM_WORLD); 
	cout <<position<<endl;
  } 
  else  /* RECEIVER CODE */ 
  {
	MPI_Recv( &position, 1 , MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	cout<<position<<endl;
	MPI_Recv( buff, position , MPI_PACKED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
	
	MPI_Finalize();
	return 0;
}
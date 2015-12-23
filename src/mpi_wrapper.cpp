#include "mpi_wrapper.h"

void MpiWorker_t::SyncAtomBool(bool& x, int root)
{
  char y;
  if(rank()==root)
	y=x;
  MPI_Bcast(&y, 1, MPI_CHAR, root, Communicator);
  x=y;
}
void MpiWorker_t::SyncVectorBool(vector< bool >& x, int root)
{
  vector <char> y;
  if(rank()==root)
	y.assign(x.begin(),x.end());
  SyncContainer(y, MPI_CHAR, root);
  if(rank()!=root)
	x.assign(y.begin(),y.end());
}
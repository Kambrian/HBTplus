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

void MpiWorker_t::SyncVectorString(vector< string >& x, int root)
{
  string buffer;
  
  if(rank()==root)
  {
    ostringstream file;
    for(auto &s: x)
      file<<s<<'\n';
    buffer=file.str();
  }
  
  SyncContainer(buffer, MPI_CHAR, root);
  
  if(rank()!=root)
  {
    x.clear();
    istringstream file(buffer);
    string s;
    while(getline(file, s))
      x.push_back(s);
  }
}

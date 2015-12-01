using namespace std;
#include <iostream>
#include <string>
#include "../boost_mpi.h"

int main(int argc, char **argv)
{
  mpi::environment env;
  mpi::communicator world;
  
  vector <vector <int> > Send(world.size()), Receive(world.size());
  for(int i=0;i<world.size();i++)
	  Send[i].assign(i, world.rank());

  
  VectorAllToAll(world, Send, Receive, MPI_INT);
  for(int i=0;i<world.size();i++)
  {
	cout<<world.rank()<<": send "<<Send[i].size()<<" elements to "<<i<<": ";
	copy(Send[i].cbegin(), Send[i].cend(), ostream_iterator<int>(cout, ", "));
	cout<<endl;
  }

  for(int i=0;i<world.size();i++)
  {
	cout<<world.rank()<<": receive "<<Receive[i].size()<<" elements from "<<i<<": ";
	copy(Receive[i].cbegin(), Receive[i].cend(), ostream_iterator<int>(cout, ", "));
	cout<<endl;
  }
  
  return 0;
}
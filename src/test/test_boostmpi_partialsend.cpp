/* it is valid to only serialize part of the struct */
using namespace std;
#include <iostream>
#include <string>
#include "mpi.h"

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;

struct TestStruct_t
{
  int a;
  vector <int> c;
  double b;
  TestStruct_t(): a(0), b(0), c(2,1)
  {
  }
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
	ar & a;
	ar & b;
// 	ar & c;
  }
};
BOOST_IS_MPI_DATATYPE(TestStruct_t)
std::ostream& operator << (std::ostream& o, vector <TestStruct_t> &x)
{
  for(auto &&a: x)
  {
   o << "[" << a.a << "; " << a.b << "; ";   
  
   auto &vec=a.c;
  copy(vec.cbegin(), vec.cend(), ostream_iterator<int>(o, ", "));
   
  o<< "]  ";
  }
  
   return o;
};
int main(int argc, char **argv)
{
  mpi::environment env;
  mpi::communicator world;
    
#define MSG_LEN 3
  vector <TestStruct_t> sendbuf(MSG_LEN), recvbuf;
  
  if(world.rank()==0)
  {
	sendbuf[1].a=1;
  sendbuf[2].b=2;
  sendbuf[1].c.push_back(-1);
  cout<<sendbuf<<endl;
	world.send(1, 0, sendbuf);
  }
  else if(world.rank()==1)
	world.recv(0, 0, recvbuf);
  
  if(world.rank()==1)
	cout<<recvbuf.size()<<" elements received: "<<recvbuf<<"\n";

  return 0;
}
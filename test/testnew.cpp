#include <iostream>
#include <new>
// #include "../datatypes.h"
class GB
{
  long data[1024*1024*1024];
};
int main()
{
  const int x=1;
  const int &y=x;
  int *p0=new int[10];
  int *p1=new (p0+5) int[3];
//   delete [] p1;
//   std::cout<<"p1 deleted\n";
  delete [] p0;
  std::cout<<"p0 deleted\n";
  std::cout<<p0<<'\n'; //the variable is still there, just content cleared.
  
  try{
  GB *p=new GB[100];
  }
  catch(std::bad_alloc &e)//this does not work. the exception is handeled before you catch it. (in C++11? should work in C++89)
  {
	std::cerr<<"failed to allocate p, "<<e.what()<<std::endl;
  }
  GB *p2=new GB[100];
  std::cout<<"success\n";
  return 0;
}
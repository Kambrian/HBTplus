//to compute the virial sizes, rmax and vmax of the host halos, using the ComovingMostBoundPosition of central subhalos as reference center
#include <cmath>
#include <iostream>
#include <string>

#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;

int main(int argc, char **argv)
{
  struct rlimit vmem, heap;
  getrlimit(RLIMIT_AS, &vmem);
  getrlimit(RLIMIT_DATA, &heap);
  cout<<"virtual limit: "<<vmem.rlim_cur<<","<<vmem.rlim_max<<endl;
  cout<<"data limit: "<<heap.rlim_cur<<", "<<heap.rlim_max<<endl<<flush;
  
  size_t s=atoll(argv[1]);
//   size_t s=322486272000LL;
  void *p=malloc(s);
  if(p==NULL)
    cerr<<"fail to allocate "<<s<<" bytes"<<endl;
  else
    cout<<"success"<<endl;
  free(p);
  
  return 0;
}

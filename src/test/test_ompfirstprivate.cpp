#include <iostream>
#include <omp.h>
#include <cstdio>
using namespace std;

class myclass
{
public:
  int x,y,z;
  myclass(): x(0),y(1) //z is uninitialized. can be anything.
  {
	printf("default called from %d\n", omp_get_thread_num());
  }
  myclass(int a, int b, int c):x(a),y(b),z(c)
  {
	printf("initialized from %d\n", omp_get_thread_num());
  }
};
int main()
{
  #pragma omp parallel num_threads(3)
  {
	int x=omp_get_thread_num();
#pragma omp for
	for(int i=0;i<3;i++)
	{
	  static myclass a(1,2,3);//static var is initialized only once!
	  printf("%d, %d: %d, %d\n", omp_get_thread_num(), i, x, a.x);
	  a.x=i;
	  printf("%d, %d: %d\n", omp_get_thread_num(), i, a.x);
	}
  }
  return 0;
}
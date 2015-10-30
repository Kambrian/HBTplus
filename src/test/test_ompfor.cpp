#include <iostream>
#include <cstdio>
#include <omp.h>
using namespace std;

int main()
{
#pragma omp parallel num_threads(3)
  for(int i=0;i<4;i++)//each thread will loop through i individually
  {
	int ithread=omp_get_thread_num();
	printf("Thread %d: i=%d\n", ithread,i);
	if(1==ithread) sleep(2);
// #pragma omp barrier
#pragma omp for //this will start to create tasks
	for(int j=0;j<i;j++){//the team of threads will loop through j together, distributing j-tasks. 
	  printf("Thread %d: i=%d, j=%d\n", ithread,i,j);
	  //there is an implicit synchronization here due to omp for.
	}
  }//in the end, each thread runs over its own i and a distributed j. If there are not enough j's, then some thread will not be assigned a j work
  return 0;
}
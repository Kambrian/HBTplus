#include <iostream>
#include <omp.h>
#include <cstdio>
using namespace std;


int main(int argc, char **argv)
{
int n=atoi(argv[1]);
//omp_set_nested(0);
omp_set_max_active_levels(1); //max_active_level 0: no para; 1: single layer;
cout<<"InPara:"<<omp_in_parallel()<<", Level:"<<omp_get_level()<<", nThreads:"<<omp_get_num_threads()<<endl;
cout<<"-----------------------------------\n";
#pragma omp parallel num_threads(2) if(n)
{
#pragma omp single
  cout<<"InPara:"<<omp_in_parallel()<<", Level:"<<omp_get_level()<<", nThreads:"<<omp_get_num_threads()<<endl;
#pragma omp parallel num_threads(3)
  {
#pragma omp single
  cout<<"InPara:"<<omp_in_parallel()<<", Level:"<<omp_get_level()<<", nThreads:"<<omp_get_num_threads()<<endl;
  }
}

  return 0;
}

#include <iostream>
#include <cstdio>
#include <omp.h>
#include <vector>
using namespace std;
int main(int argc, char **argv)
{
  int n=3;
  if(argc>1) n=atoi(argv[1]);
#pragma omp parallel num_threads(3)
  {
	static vector <int> x(3, 1); //by default, x is private, and is initialized properly
	int y[n];//this is allowed by C99 and g++
	for(int i=0;i<n;i++) y[i]=1;
	x[0]=-omp_get_thread_num();
	y[0]=x[0];
#pragma omp barrier
#pragma omp critical
	{
	  cout<<"("<<x[0]<<","<<x[1]<<","<<x[2]<<")  ";
	  cout<<"("<<y[0]<<","<<y[1]<<","<<y[2]<<")"<<endl;
	}
  }
  return 0;
}
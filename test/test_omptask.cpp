#include <iostream>
#include <cstdio>
#include <omp.h>
#include <vector>
#include <ctime>
#define N 10000
using namespace std;
int work(int i, int j)
{
  int x=i+j*2+j*j*j+i*i;
  x*=x;
  return x;
}
int main(int argc, char **argv)
{
  time_t t0,t1;
  
      t0=time(NULL);
#pragma omp parallel num_threads(4)
  {
#pragma omp for
	for(int i=0;i<N;i++)
	  for(int j=0;j<i;j++)
		int x=work(i,j);
  }
  t1=time(NULL);
  cout<<"outer for:"<<t1-t0<<endl; 
  
  t0=time(NULL);
#pragma omp parallel num_threads(4)
  {
	for(int i=0;i<N;i++)
#pragma omp for
	  for(int j=0;j<i;j++)
		int x=work(i,j);
  }
  t1=time(NULL);
  cout<<"inner for:"<<t1-t0<<endl;
  
      t0=time(NULL);
	  int x;
	for(int i=0;i<N;i++)
	  for(int j=0;j<i;j++)
		x+=work(i,j);
  t1=time(NULL);
  cout<<"serial:"<<t1-t0<<endl;
  
  
  t0=time(NULL);
#pragma omp parallel num_threads(4)
  {
#pragma omp single nowait
	for(int i=0;i<N;i++)
	  for(int j=0;j<i;j++)
#pragma omp task firstprivate(i,j)
		int x=work(i,j);
#pragma omp taskwait		
  }
  t1=time(NULL);
  cout<<"task nowait: "<<t1-t0<<endl;//slow, because it's too fine-grained..

#pragma omp parallel num_threads(4)
  {
#pragma omp single
	for(int i=0;i<N;i++)
	  for(int j=0;j<i;j++)
#pragma omp task firstprivate(i,j)
		int x=work(i,j);
  }
  t1=time(NULL);
  cout<<"task: "<<t1-t0<<endl;//slowest!
  
  return 0;
}
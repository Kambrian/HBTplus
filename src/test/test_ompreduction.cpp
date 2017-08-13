#include <iostream>
#include <cstdio>
#include <omp.h>
using namespace std;

#define N 100
int main()
{
  int x[N];
  for(int i=0;i<N;i++)
    x[i]=i;
  swap(x[0], x[30]);
  int imin,imax;
  imin=imax=x[0];
#pragma omp parallel for reduction(max: imax), reduction(min: imin)
  for(int i=1;i<N;i++)
  {
    if(x[i]>imax)
      imax=x[i];
    if(x[i]<imin) //cannot use elseif here, because of the initialization of reduction vars
      imin=x[i];
  }
  
  cout<<imin<<","<<imax<<endl;
  return 0;
}
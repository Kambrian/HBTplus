#include <iostream>
#include <cstdio>
#include <omp.h>
#include <random>
#include <sstream>

using namespace std;

int main()
{
  int N=500;
  float x[2][N];
#pragma omp parallel num_threads(2) 
  {
    int ithread=omp_get_thread_num();
    default_random_engine generator(ithread+100);
    uniform_real_distribution<float> distribution(0.0,1.0);
    stringstream ss;
    ss<<ithread<<":"<<distribution(generator)<<endl;
    cout<<ss.str();
    for(int i=0;i<N;i++)//each thread will loop through i individually
    {
      x[ithread][i]=distribution(generator);
    }
  }
  for(int i=0;i<2;i++)
  {
    cout<<"[ ";
    for(int j=0;j<N;j++)
      cout<<x[i][j]<<",";
    cout<<"]"<<endl;
  }
  return 0;
}
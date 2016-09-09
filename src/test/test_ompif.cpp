#include <iostream>
#include <cstdio>
#include <omp.h>
#include <chrono>
#include <vector>
using namespace std;

class Timer_t
{
public:
  vector <chrono::high_resolution_clock::time_point> tickers;
  Timer_t()
  {
	tickers.reserve(20);
  }
  void Tick()
  {
	tickers.push_back(chrono::high_resolution_clock::now());
  }
  void Reset()
  {
	tickers.clear();
  }
  size_t Size()
  {
	return tickers.size();
  }
  double GetSeconds(int itick)
  /*get the time spent from the previous tick to the current tick*/
  {
	return GetSeconds(itick, itick-1);
  }
  double GetSeconds(int itick, int itick0)
  /*get the time spent from itick0 to itick*/
  {
	double tsign=1.;
	if(itick<itick0)
	{
	  tsign=-1.;
	  swap(itick, itick0);
	}
	if(itick>Size()) itick=Size();
	if(itick0<0) itick0=0;
	return tsign*chrono::duration_cast<chrono::duration<double> >(tickers[itick]-tickers[itick0]).count();
  }
};

int Acker(int m, int n)
{
  if(m==0) return n+1;
  if(m>0&&n==0) return Acker(m-1, 1);
  if(m>0&&n>0) return Acker(m-1, Acker(m, n-1));
  return 0;
}
void work()
{
  int x=Acker(4,1);
}
int main(int argc, char **argv)
{
  int flag=atoi(argv[1]);
  int n=atoi(argv[2]);
  Timer_t timer;
  timer.Tick();
#pragma omp parallel for if(flag)
  for(int i=0;i<n;i++)
  {
    work();
  }
  timer.Tick();
#pragma omp parallel for
  for(int i=0;i<n;i++)
    work();
  timer.Tick();
  cout<<"Time: "<<timer.GetSeconds(1)<<","<<timer.GetSeconds(2)<<endl;
  return 0;
}
#include <iostream>
#include <cstdio>
#include <omp.h>
using namespace std;

class A
{
  int i, j;
  void task1()
  {
    printf("1, %d, %d\n", i, j);
  }
  void task2()
  {
    printf("2, %d, %d\n", i, j);
  }
  void task3()
  {
    printf("3, %d, %d\n", i, j);
  }
public:
  A(int i, int j): i(i), j(j)
  {
  }
  void job()
  {
#pragma omp parallel
#pragma omp single
    {
#pragma omp task
      task1();
#pragma omp task
      task2();
#pragma omp task
      task3();
    }
  }

};

int main(int argc, char **argv)
{
  A a(10, 20);
  a.job();
  
  return 0;
}
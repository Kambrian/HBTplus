#include <iostream>
#include <new>
#include <vector>
// #include "../datatypes.h"

int main()
{
  std::vector <int> x;
  int i;
  i=x.size()-1;
  std::cout<<i<<","<<x.size()-1<<std::endl;
  //output: -1,18446744073709551615 (int and size_type)
  return 0;
}

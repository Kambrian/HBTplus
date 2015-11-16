#include <iostream>
#include <vector>
#include "../mymath.h"

int main(int argc, char **argv)
{
#define test(x,y)  cout<<x<<","<<y<<" : "<<ClosestFactors(x,y)<<endl;;
  /*
test(10,3);
test(5,2);
test(8,3);
test(20,1);
test(20,2);
*/
  
  test(atoi(argv[1]), atoi(argv[2]));
  
  return 0;
}
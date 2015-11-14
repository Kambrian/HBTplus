#include <iostream>
#include <vector>
#include <numeric>
using namespace std;
int main()
{
  vector <int> x(10,1);
  int s=0;
  s=accumulate(x.begin(),x.end(),s);
  cout<<s<<endl;
  
  s=accumulate(x.begin()+1,x.begin()+5, s);
  cout<<s<<endl;
  
  
  int init = 100;
  int numbers[] = {10,20,30};

  std::cout << "using default accumulate: ";
  std::cout << std::accumulate(numbers,numbers+3,init);
  std::cout << '\n';
  
  return 0;
}
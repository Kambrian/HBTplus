#include <iostream>
#include <vector>

using namespace std;
int main()
{
  vector <int> x(3,0);
  for(auto &&a: x)
    a=10;
  for(auto &a:x)
    cout<<a<<",";
  cout<<endl;
  
  for(int i=1;i<1;i++)
    cout<<"i="<<i<<endl;
  return 0;
}
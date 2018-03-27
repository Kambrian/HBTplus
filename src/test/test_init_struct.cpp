#include <iostream>
#include <vector>
using namespace std;

struct myclass
{
  int x,y,z;
  int a[3];
  void print()
  {
    cout<<x<<','<<y<<','<<z<<','<<a[0]<<','<<a[1]<<','<<a[2]<<endl;
  }
};
int main()
{
vector <myclass> arr(10);
myclass x;
x.print();//everything un-initialized
for(auto &&p: arr) //the vector elements appear to be all value initialized
  p.print();

  return 0;
}
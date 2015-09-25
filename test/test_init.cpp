#include <iostream>
using namespace std;

class myclass
{
public:
  int x,y,z;
  myclass(): x(0),y(1) //z is uninitialized. can be anything.
  {
  }
};
int main()
{
  myclass c=myclass(),d;
  cout<<c.x<<','<<c.y<<','<<c.z<<endl;
  cout<<d.x<<','<<d.y<<','<<d.z<<endl;
  /*output:
    0,1,4196688
	0,1,4196128
	*/
  return 0;
}
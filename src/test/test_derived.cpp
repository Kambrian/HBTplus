#include <iostream>
#include <cstdio>
using namespace std;

class myclass
{
public:
 int x;
  myclass(): x(0)
  {
	printf("default called from base class\n");
  }
  virtual void myprint()
  {
    x=1;
    printf("myprint from base class,x=%d\n",x);
  }
  void print()
  {
    myprint();
  }
};
class yourclass: public myclass
{
  void myprint()
  {
    x=2;
    printf("myprint from drived class, x=%d\n",x);
  }
};

int main()
{
myclass a;
yourclass b;

a.print();
b.print();

return 0;
}
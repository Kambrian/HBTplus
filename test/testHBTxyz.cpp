#include <iostream>
#include <array>
// #include "../datatypes.h"
template <class T>
class XYZ
{
  std::array<T,3> data;
public:
  XYZ(T x=0, T y=0, T z=0) //:data{x,y,z}
  {
	data[0]=x;
	data[1]=y;
	data[2]=z;
  }
  XYZ(const XYZ &s)
  {
	std::cout<<"assign"<<std::endl;
	data=s.data;
  }
  T & operator [](int i)
  {
	return data[i];
  }
  ~XYZ()
  {
	std::cout<<"deleted"<<std::endl;
  }
};
template <class T>
std::ostream& operator << (std::ostream& o, XYZ <T> &a)
{
  
   o << "(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
   return o;
};
typedef XYZ <float> HBTxyz;
class myxyz: public XYZ <double>
{
public:
  myxyz(double x,double y, double z): XYZ <double> (x,y,z)
  {
  }
  //default destructor is auto called.
};
int main()
{
  XYZ <float> x={1.,2.,3.};
  std::cout<<x[0]<<","<<x[1]<<","<<x[2]<<std::endl;
  XYZ <float> y(x);
  HBTxyz &z=y;
  myxyz b(1.,2.,3.);
  std::cout<<y<<std::endl;
  y[0]=0.;
  std::cout<<x<<y<<z<<b<<"\n";
  std::cout<<"loop test:\n";
  for(int i=0;i<5;i++)
  {
	XYZ <float> x;
	std::cout<<i<<": ";
  }
  return 0;
}
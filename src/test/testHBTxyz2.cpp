#include <iostream>
#include <string.h>
#include <vector>
#include "../datatypes.h"
template <class T>
class XYZ
{
  T data[3];
public:
  XYZ(T x, T y, T z) //:data{x,y,z}
  {
	data[0]=x;
	data[1]=y;
	data[2]=z;
  }
  XYZ(const XYZ &s)
  {
	std::cout<<"assign"<<std::endl;
	memcpy(data, s.data, sizeof(T)*3);
  }
  T * address()
  {
	return data;
  }
  T & operator [](int i)
  {
	return data[i];
  }
  ~XYZ()
  {
// 	delete data;
	std::cout<<"deleted"<<std::endl;
  }
};
template <class T>
std::ostream& operator << (std::ostream& o, XYZ <T> &a)
{
  
   o << "(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
   return o;
};
std::ostream& operator << (std::ostream& o, HBTxyz &a)
{
  
   o << "(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
   return o;
};
typedef double dxyz[3];
std::ostream& operator << (std::ostream& o, dxyz &a)
{
  
   o << "(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
   return o;
};

int main()
{
  XYZ <float> x(1.,2.,3.);
  std::cout<<x[0]<<","<<x[1]<<","<<x[2]<<std::endl;
  XYZ <float> y(x);
  XYZ <float> &z=y;
  HBTxyz a={-1,-2,-3}, b={-4,-5,-6};
  std::cout<<y<<z<<std::endl;
  y[0]=0.;
  std::cout<<x<<y<<z<<a<<"\n";
  
  a=b;
  cout<<a<<b<<&a<<","<<&b<<endl;
  vector <HBTxyz> v;
  float data[]={1,2,3,4,5,6};
//   v.assign((float (*)[3])data, (float (*)[3])data+2);
  v.resize(2);
  copy(data, data+6, (HBTReal *) v.data());
  cout<<v[0]<<v[1]<<endl;
  cout<<&a<<','<<a.data()<<endl;
  return 0;
}
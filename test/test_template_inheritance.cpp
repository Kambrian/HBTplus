#include <iostream>
#include <cstring>

typedef int HBTInt;

template <class T>
class List_t
{ 
  HBTInt N;
public:
  T * Data;
  List_t(HBTInt n=0, T *data=NULL)
  {
	N=n;
	if(data)//memory can be shared, not always allocated. ToDo: implement shared_ptr?
	  Data=data;
	else
	  Data=new T[n];
  }
  T & operator [](HBTInt i)
  {
	return Data[i];
  }
  HBTInt size()
  {
	return N;
  }
  void clear()//the user is responsible for cleaning up.
  {
	delete [] Data;
	N=0;
  }
};
struct ParticleReference_t
{
  HBTInt Id;
  HBTInt Index;
};
template <class T>
class ParticleList_t: public List_t <T>
{
public:
  T & operator [](ParticleReference_t ref) //overload for reference
  {
	return List_t <T>::Data[ref.Index]; //or this->Data could also work; but Data along is unqualified.
  }
};
int main()
{
  List_t <ParticleReference_t> L(10);
  std::cout<<L.size()<<std::endl;
}


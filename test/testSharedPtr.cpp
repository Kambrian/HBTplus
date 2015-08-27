#include <memory>

template <class T>
class SharedList_t
{
  HBTInt N;
  std::shared_ptr <T> Data;
public:
  SharedList_t(): N(0), Data(NULL)
  {
  }
  SharedList_t(HBTInt n)
  {
	Data=new T[n];
	N=n;
  }
  SharedList_t(HBTInt n, void *data, std::shared_ptr <T> storage):
  N(n), Data(storage, data)
  {
  }
  //the default copy constructor would be shallow copy.
  SharedList_t(SharedList_t & list): N(list.N), Data(list.Data)
  {
  }
  void deep_copy(SharedList_t & list)
  {//deep copy
	N=list.N;
	Data=new T[N];
	memcpy(Data, list.Data, sizeof(T)*N);
  }
  void reset()
  {
	N=0;
	Data.reset();
  }
  HBTInt size()
  {
	return N;
  }
  T & operator [] (HBTInt i)
  {
	return ((T *)Data)[i];
  }
};
struct HaloCatalogue
{  
  SharedList_t <int> * Halos;
  int Nhalos;
}
int main()
{
  SharedList_t <int> H1(10),H2(5,&(H1[5]),&H1);
  return 0;
}
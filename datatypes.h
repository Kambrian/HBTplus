#ifndef DATATYPES_INCLUDED

#include <iostream>
#include <cstring>
// #include <memory>

/*datatype for input particle data*/
#ifdef INPUT_REAL8
typedef double IDatReal;
#else
typedef float IDatReal;
#endif

/*datatype for input particle IDs*/
#ifdef INPUT_INT8
typedef long IDatInt;
#else
#ifdef INPUT_UINT4
typedef unsigned IDatInt;
#else
typedef int IDatInt;
#endif
#endif

/*datatype for internal calculation and output*/
#ifdef HBT_REAL8
typedef double HBTReal;
#else
typedef float HBTReal;
#endif

// the user should ganrantee that HBTInt can at least hold NP_DM (when HBTPID_RANKSTYLE is defined)
#ifdef HBT_INT8
typedef long HBTInt;  
#define HBTIFMT "%ld"
#else 
typedef int HBTInt;
#define HBTIFMT "%d"
#endif

#if (defined INPUT_REAL8 && defined HBT_REAL8)||(!defined INPUT_REAL8 && !defined HBT_REAL8)
 #define SAME_REALTYPE
#endif

#if (defined INPUT_INT8 && defined HBT_INT8)||(!defined INPUT_INT8 && !defined HBT_INT8)
 #define SAME_INTTYPE
#endif

//constants
#define FRSH_GRPCAT -1
#define FRSH_SUBCAT -2
#define FRSH_SRCCAT -3
//#define FRSH_MBDCAT -4

#define GROUP_FORMAT_GADGET4 40
#define GROUP_FORMAT_GADGET3_INT 30
#define GROUP_FORMAT_GADGET3_LONG 31

typedef HBTReal HBTxyz[3];  //3-d pos/vel data

namespace SpecialConst
{
  const HBTInt NullParticleId=-100;
  const HBTInt NullSnapshotId=-100;
  const HBTInt NullHaloId=-100;
  const HBTInt NullSubhaloId=-100;
  const HBTInt NullTrackId=-100;
  
  const HBTInt UniverseHaloId=-200; //the universe as a halo, to contain particles that does not belong to any halo
  
  const HBTxyz NullCoordinate={0.,0.,0.};
//   const Particle_t NullParticle(NullParticleId, NullParticleId, NullCoordinate, NullCoordinate);
};

struct Particle_t
{
  HBTInt ParticleId;
  HBTInt ParticleIndex;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
};
//TODO: differentiate between clean and dirty list, or initiate and follow list.
template <class T>
class List_t
/* a base class for list; only empty constructor*/
{ 
protected:  
  HBTInt N;
  T * Data; //this is only copied. never allocated by itself.
public:
  List_t(): N(0), Data(nullptr)
  {
  }
  T & operator [](const HBTInt index)
  {
	return Data[index];
  }
  HBTInt Size()
  {
	return N;
  }
};
template <class T>
class ShallowList_t: public List_t <T>
/* a list whose memory is only a link pointing to already allocated buffer.
 * It involves no new or delete inside its constructor or destructor.
 * so it is only a "virtual" list.*/
{
  typedef List_t <T> BaseList_t;
public:
  ShallowList_t(): BaseList_t()
  {
  }
  ShallowList_t(const HBTInt n, T * const data): BaseList_t::N(n), BaseList_t::Data(data)
  {
  }
  ShallowList_t(const ShallowList_t &l): BaseList_t::N(l.N), BaseList_t::Data(l.Data)
  {
  }
  void Bind(HBTInt n, T *data)
  {
	BaseList_t::N=n;
	BaseList_t::Data=data;
  }
};
template <class T>
class DeepList_t: public List_t <T>
/* DeepList_t manages memory actively. automatically allocate memory upon construction;
 automatically frees memory upon destruction;*/
{
typedef List_t <T> BaseList_t;  
public:
  DeepList_t(): BaseList_t ()
  {
  }
  template <class U>
  DeepList_t(const DeepList_t<U> & l)
  {
	Fill(l.n, l.Data);
  }
  void Resize(const HBTInt n)
  {
	if(n!=BaseList_t::N)
	{
	  delete [] BaseList_t::Data;
	  BaseList_t::N=n;
	  BaseList_t::Data=new T[n];
	}
  }
  template <class U>
  void Fill(const HBTInt n, U * const data)
  {
	Resize(n);
	for(HBTInt i=0;i<n;i++)
	  BaseList_t::Data[i]=data[i];
  }
  void Bind(const HBTInt n, T * const data)
  /*bind data to DeepList's internal pointer, instead of copying the data.
   *Caution: data has to be created with new T[]; otherwise may cause segfault upon destruction. use with care. */
  {
	delete [] BaseList_t::Data;
	BaseList_t::Data=data;
	BaseList_t::N=n;
  }
  void Clear()
  {
	delete [] BaseList_t::Data;
	BaseList_t::Data=nullptr;
	BaseList_t::N=0;
  }
  ~DeepList_t()
  {
	delete [] BaseList_t::Data;
// 	std::cout<<"List cleaned\n";
// 	BaseList_t::N=0; //no need to do this. the object is out of scope (destroyed) after this is called.
  }
};
#define DATATYPES_INCLUDED
#endif

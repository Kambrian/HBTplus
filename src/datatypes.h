#ifndef DATATYPES_INCLUDED

#include <iostream>
#include <cstring>
using namespace std;
#include <array>
#include <vector>
// #include <memory>
#ifdef DM_ONLY
#undef UNBIND_WITH_THERMAL_ENERGY
#undef HAS_THERMAL_ENERGY
#endif

#ifdef UNBIND_WITH_THERMAL_ENERGY
#define HAS_THERMAL_ENERGY
#endif

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

// typedef HBTReal HBTxyz[3];  //3-d pos/vel data
/*inline void copyHBTxyz(HBTxyz & dest, const HBTxyz & src)
{
  memcpy(dest, src, sizeof(HBTxyz));
}*/
typedef array <HBTReal, 3> HBTxyz;
inline void copyHBTxyz(HBTxyz &dest, const HBTxyz &src)
{
  /*copy for std:arr implementation*/
  dest=src;
}
template <class T>
inline void copyHBTxyz(HBTxyz &dest, const T src[3])
{
  dest[0]=src[0];
  dest[1]=src[1];
  dest[2]=src[2];
}
template <class T1, class T2>
inline void copyXYZ(T1 & dest, const T2 src)
{
  dest[0]=src[0];
  dest[1]=src[1];
  dest[2]=src[2];
}

namespace SpecialConst
{
  const HBTInt NullParticleId=-1;
  const HBTInt NullSnapshotId=-1;
  const HBTInt NullHaloId=-1;//do not change this.
  const HBTInt NullSubhaloId=-1;
  const HBTInt NullTrackId=-1;
    
  const HBTxyz NullCoordinate={0.,0.,0.};
//   const Particle_t NullParticle(NullParticleId, NullParticleId, NullCoordinate, NullCoordinate);
};

template <class T>
class VectorView_t
/* similar to vector, but never actively manage memory; only bind to existing memory*/
{
public:
  HBTInt N;
  T * Data; //this is only copied. never allocated by itself.
  VectorView_t(): N(0), Data(nullptr)
  {
  }
  VectorView_t(const HBTInt n, T * const data): N(n), Data(data)
  {
  }
  void Bind(const HBTInt n, T * const data)
  {
	N=n;
	Data=data;
  }
  void Bind(T * const data)
  {
	Data=data;
  }
  void ReBind(const HBTInt n)
  {
	N=n;
  }
  void IncrementBind()
  {
	N++; 
  }
  T * data() const
  {
	return Data;
  }
  T & operator [](const HBTInt index) const
  {
	return Data[index];
  }
  HBTInt size() const
  {
	return N;
  }
  void PushBack(T x)
  /*memory is never reallocated*/
  {
	Data[N]=x;
	N++;
  }
  T * begin()
  {
	return Data;
  }
  T* end()
  {
	return Data+N;
  }
  T & back()
  {
	return Data[N-1];
  }
};

typedef enum
{
  TypeGas=0,
  TypeDM,
  TypeDisk,
  TypeBulge	,
  TypeStar,
  TypeBndry,
  TypeMax
} ParticleType_t;

struct LocatedParticle_t
{
  HBTInt index;
  HBTReal d2; //distance**2
  LocatedParticle_t(){};
  LocatedParticle_t(HBTInt index, HBTReal d2):index(index),d2(d2)
  {}
};
inline bool CompLocatedDistance(const LocatedParticle_t &a, const LocatedParticle_t &b)
{
  return a.d2<b.d2;
}
class ParticleCollector_t
{
public:
  virtual void Collect(HBTInt index, HBTReal d2)=0;
};
class LocatedParticleCollector_t: public ParticleCollector_t
//a simple collector that appends all found particles to a vector.
{
public:
  vector <LocatedParticle_t> Founds;
  LocatedParticleCollector_t(HBTInt n_reserve=0): Founds()
  {
    Founds.reserve(n_reserve);
  }
  void Collect(HBTInt index, HBTReal d2)
  {
    Founds.emplace_back(index, d2);
  }
  void Clear()
  {
    Founds.clear();
  }
};
class NearestNeighbourCollector_t: public ParticleCollector_t
//a collector for nearest neighbour search. keeps the nearest neighbour particle.
{
public:
  HBTInt Index;
  HBTReal D2;//distance squared
  NearestNeighbourCollector_t():Index(0), D2(-1)
  {}
  void Collect(HBTInt index, HBTReal d2)
  {
    if(d2<D2)
    {
      D2=d2;
      Index=index;
    }
  }
  bool IsEmpty()
  {
    return D2<0.;
  }
};
#define DATATYPES_INCLUDED
#endif

#ifndef DATATYPES_INCLUDED

#include <iostream>
#include <iterator>
#include <cstring>
using namespace std;
#include <array>
// #include <memory>
#ifdef DM_ONLY
#undef UNBIND_WITH_THERMAL_ENERGY
#undef HAS_THERMAL_ENERGY
#endif

#ifdef UNBIND_WITH_THERMAL_ENERGY
#ifndef HAS_THERMAL_ENERGY
#define HAS_THERMAL_ENERGY
#endif
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
#define MPI_HBT_REAL MPI_DOUBLE
#else
typedef float HBTReal;
#define MPI_HBT_REAL MPI_FLOAT
#endif

// the user should ganrantee that HBTInt can at least hold NP_DM 
#ifdef HBT_INT8
typedef long HBTInt;  
#define HBTIFMT "%ld"
#define MPI_HBT_INT MPI_LONG
#else 
typedef int HBTInt;
#define HBTIFMT "%d"
#define MPI_HBT_INT MPI_INT
#endif

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
  const HBTInt NullParticleId=-1;//reserved special id, should not be used by input simulation data
  const HBTInt NullSnapshotId=-1;
  const HBTInt NullHaloId=-1;//do not change this.
  const HBTInt NullSubhaloId=-1;
  const HBTInt NullTrackId=-1;
    
  const HBTxyz NullCoordinate={0.,0.,0.};
//   const Particle_t NullParticle(NullParticleId, NullParticleId, NullCoordinate, NullCoordinate);
};

struct IdRank_t
{
  HBTInt Id;
  int Rank;
  IdRank_t(){};
  IdRank_t(HBTInt id, int rank): Id(id), Rank(rank)
  {
  }
};
#ifdef HBT_INT8
#define MPI_HBTRankPair MPI_LONG_INT
#else
#define MPI_HBTRankPair MPI_2INT
#endif
inline bool CompareRank(const IdRank_t &a, const IdRank_t &b)
{
  return (a.Rank<b.Rank);
}

template <class T>
class VectorView_t
/* similar to vector, but never actively manage memory; only bind to existing memory*/
{
public:
  typedef T * iterator;
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

enum ParticleType_t:int
{
  TypeGas=0,
  TypeDM,
  TypeDisk,
  TypeBulge	,
  TypeStar,
  TypeBndry,
  TypeMax
};

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

#define DATATYPES_INCLUDED
#endif

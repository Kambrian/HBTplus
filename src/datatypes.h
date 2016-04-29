#ifndef DATATYPES_INCLUDED

#include <iostream>
#include <iterator>
#include <cstring>
using namespace std;
#include <array>
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
#define MPI_HBT_REAL MPI_DOUBLE
#else
typedef float HBTReal;
#define MPI_HBT_REAL MPI_FLOAT
#endif

// the user should ganrantee that HBTInt can at least hold NP_DM (when HBTPID_RANKSTYLE is defined)
#ifdef HBT_INT8
typedef long HBTInt;  
#define HBTIFMT "%ld"
#define MPI_HBT_INT MPI_LONG
#else 
typedef int HBTInt;
#define HBTIFMT "%d"
#define MPI_HBT_INT MPI_INT
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
#define GROUP_FORMAT_GADGET2_INT 20
#define GROUP_FORMAT_GADGET2_LONG 21

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

struct IdRank_t
{
  HBTInt Id;
  int Rank;
  IdRank_t()=default;
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
#define DATATYPES_INCLUDED
#endif

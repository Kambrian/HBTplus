/* CAUTION: if your data exceeds SIZET_MAX, there would be problem allocating memory*/
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
#define HBTIFMT "%lld"
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

#define FLAG_LOAD_ID  0b001
#define FLAG_LOAD_POS 0b010
#define FLAG_LOAD_VEL 0b100

//auxiliary datatype
template <class T>
class XYZ
{
  T data[3];
public:
  XYZ(T x=0., T y=0., T z=0.) //:data{x,y,z}
  {
	data[0]=x;
	data[1]=y;
	data[2]=z;
  }
  XYZ(const XYZ <T> &s)
  {
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
};
template <class T>
std::ostream& operator << (std::ostream& o, XYZ <T> &a)
{
   o << "(" << a[0] << ", " << a[1] << ", " << a[2] << ")";
   return o;
};
typedef XYZ <HBTReal> HBTxyz;
typedef XYZ <float> Float_xyz;
typedef XYZ <double> Double_xyz;
// typedef HBTReal HBTxyz[3];  //3-d pos/vel data

namespace SpecialConst
{
  const HBTInt NullParticleId=-100;
  const HBTInt NullSnapshotId=-100;
  const HBTInt NullHaloId=-100;
  const HBTInt NullSubhaloId=-100;
  const HBTInt NullTrackId=-100;
  
  const HBTInt UniverseHaloId=-200; //the universe as a halo, to contain particles that does not belong to any halo
  
  const HBTxyz NullCoordinate(0.,0.,0.);
//   const Particle_t NullParticle(NullParticleId, NullParticleId, NullCoordinate, NullCoordinate);
};

struct ParticleReference_t
{
  union 
  {
	HBTInt Id;
	HBTInt Index;
  };
};

class Particle_t: public ParticleReference_t
{
public:
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
};

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
  T & operator [](ParticleReference_t ref)
  {
	return Data[ref.Index];
  }
  HBTInt Size()
  {
	return N;
  }
  void Clear()//the user is responsible for cleaning up.
  {
	delete [] Data;
	N=0;
  }
};

#ifndef IniSnap //initial snapshot to start processing
#define IniSnap 0
#endif

#define DATATYPES_INCLUDED
#endif

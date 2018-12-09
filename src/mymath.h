#ifndef MYMATH_HEADER_INCLUDED
#define MYMATH_HEADER_INCLUDED

#include <iostream>
#include <iterator>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include <chrono>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <mpi.h>

#include "datatypes.h"

#define VecDot(x,y) ((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2])
#define VecNorm(x) VecDot(x,x)

extern int GetGrid(HBTReal x, HBTReal step, int dim);
extern int AssignCell(const HBTxyz & Pos, const HBTxyz &step, const vector <int> &dims);

template <class T>
void VectorFree(vector <T> &x)
{
  vector <T>().swap(x);
}

template <class T, class T2>
size_t CompileOffsets(const vector <T> &Counts, vector <T2> &Offsets)
{
  size_t offset=0;
  Offsets.resize(Counts.size());
  for(size_t i=0;i<Counts.size();i++)
  {
	Offsets[i]=offset;
	offset+=Counts[i];
  }
  return offset;
}

template <class CountIterator_t, class OffsetIterator_t>
size_t CompileOffsets(CountIterator_t CountBegin, CountIterator_t CountEnd, OffsetIterator_t OffsetBegin)
{
  size_t offset=0;
  auto it_off=OffsetBegin;
  for(auto it=CountBegin;it!=CountEnd;++it)
  {
	*it_off++=offset;
	offset+=*it;
  }
  return offset;
}

template <class T, class UnaryPredicate>
inline void RemoveFromVector(vector <T> &v, UnaryPredicate p)
{
  v.erase(remove_if(v.begin(), v.end(), p), v.end());
}

class Timer_t
{
public:
  vector <chrono::high_resolution_clock::time_point> tickers;
  Timer_t()
  {
	tickers.reserve(20);
  }
  void Tick()
  {
	tickers.push_back(chrono::high_resolution_clock::now());
  }
  void Tick(MPI_Comm comm)
  //synchronized tick. wait for all processes to tick together.
  {
	MPI_Barrier(comm);
	Tick();
  }
  void Reset()
  {
	tickers.clear();
  }
  size_t Size()
  {
	return tickers.size();
  }
  int FixTickNum(int itick)
  {
    return itick<0?itick+Size():itick;
  }
  double GetSeconds(int itick=-1)
  /*get the time spent from the previous tick to the current tick
   * if itick not specified, return the current interval
   * if itick<0, it will be interpreted as end()+itick */
  {
	itick=FixTickNum(itick);
	return GetSeconds(itick, itick-1);
  }
  double GetSeconds(int itick, int itick0)
  /*get the time spent from itick0 to itick*/
  {
    itick=FixTickNum(itick);
    itick0=FixTickNum(itick0);
    if(itick<itick0)
      swap(itick, itick0);
	
    return chrono::duration_cast<chrono::duration<double> >(tickers[itick]-tickers[itick0]).count();
  }
};

#define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){fprintf(stderr,"Error opening file %s\n",filename);	fflush(stderr); exit(1);}
// #ifdef PERIODIC_BDR
// #define NEAREST(x) (((x)>BOXHALF)?((x)-BOXSIZE):(((x)<-BOXHALF)?((x)+BOXSIZE):(x)))
			/*this macro can well manipulate boundary condition because 
			* usually a halo is much smaller than boxhalf
			* so that any distance within the halo should be smaller than boxhalf */
// #endif
#define get_bit(x,k) (((x)&(1<<k))>>k)
extern int count_pattern_files(char *filename_pattern);
// extern std::ostream& operator << (std::ostream& o, HBTxyz &a);
extern void swap_Nbyte(void *data2swap,size_t nel,size_t mbyte);
extern size_t SkipFortranBlock(FILE *fp, bool NeedByteSwap);
template <class T, std::size_t N>
ostream& operator<<(ostream& o, const array<T, N>& arr)
{
  o<<"(";
  copy(arr.cbegin(), arr.cend(), ostream_iterator<T>(o, ", "));
  o<<")";
  return o;
}

template <class T>
ostream& operator<<(ostream& o, const vector<T>& vec)
{
  o<<"(";
  copy(vec.cbegin(), vec.cend(), ostream_iterator<T>(o, ", "));
  o<<")";
  return o;
}

inline size_t fread_swap(void *buf,const size_t member_size, const size_t member_count,FILE *fp, const bool FlagByteSwap)
{
	size_t Nread;
	Nread=std::fread(buf,member_size,member_count,fp);
	if(FlagByteSwap)
	swap_Nbyte(buf,member_count,member_size);
	return Nread;
}
inline bool file_exist(char * filename)
{ struct stat buffer;   
  return (stat(filename, &buffer) == 0); 
}
inline long int BytesToEOF(FILE *fp)
{
  fpos_t fpos;
  fgetpos (fp,&fpos);
  
  long int offset=ftell(fp);
  fseek(fp, 0L, SEEK_END);
  long int offset_end=ftell (fp);
  
  fsetpos(fp, &fpos);
  
  return (offset_end-offset);
}
inline HBTReal position_modulus(HBTReal x, HBTReal boxsize)
{//shift the positions to within [0,boxsize)
	HBTReal y;
	if(x>=0&&x<boxsize) return x;
	y=x/boxsize;
	return (y-floor(y))*boxsize;
}
inline HBTReal Distance2(const HBTReal x[3], const HBTReal y[3])
{
	HBTxyz dx;
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];
	dx[2]=x[2]-y[2];
	return dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
}
inline HBTReal Distance(const HBTReal x[3], const HBTReal y[3])
{
	return sqrt(Distance2(x,y));
}
inline HBTReal Distance(const HBTxyz &x, const HBTxyz &y)
{
  return Distance(x.data(), y.data());
}
inline HBTReal Distance(const HBTxyz &x, const HBTReal y[3])
{
  return Distance(x.data(), y);
}
inline HBTReal Distance(const HBTReal x[3], const HBTxyz &y)
{
  return Distance(x, y.data());
}

template <class T>
class FortranBlock
{
  vector <T> Data;
  typedef T Txyz[3];
public:
  FortranBlock(): Data()
  {
  }
  FortranBlock(FILE *fp, const size_t n_read, const size_t n_skip, bool NeedByteSwap=false): Data(n_read)
  {
    Read(fp, n_read, n_skip, NeedByteSwap);
  }
  void Read(FILE *fp, const size_t n_read, const size_t n_skip, bool NeedByteSwap=false)
/*read n_read members from the current block of fp. 
 * skip n_skip elements before reading. 
 * T specify the input datatype. if T and U has the same size, read directly into outbuffer; otherwise the elements are converted from type U to type T in a temporary buffer and then copied to outbuffer.
 */
  {
  #define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
  #define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)
	  int blocksize,blocksize2;
	  ReadBlockSize(blocksize);
	  size_t block_member_size=sizeof(T);
	  Data.resize(n_read);
	  fseek(fp, n_skip*block_member_size, SEEK_CUR);
	  myfread(Data.data(), block_member_size, n_read, fp);
	  fseek(fp, blocksize-(n_skip+n_read)*block_member_size, SEEK_CUR);
	  ReadBlockSize(blocksize2);
	  assert(blocksize==blocksize2);
  #undef ReadBlockSize
  #undef myfread
  }
  const T * data()
  {
	return Data.data();
  }
  const T & operator [](const size_t index)
  {
	return Data[index];
  }
  HBTInt size() const
  {
	return Data.size();
  }
  T * begin()
  {
	return Data.data();
  }
  T* end()
  {
	return Data.data()+Data.size();
  }
  Txyz * data_reshape()
  {
	return (Txyz *)Data.data();
  }
};

extern int LargestRootFactor(int N, int dim);
extern vector <int> ClosestFactors(int N, int dim);
extern void AssignTasks(HBTInt worker_id, HBTInt nworkers, HBTInt ntasks, HBTInt &task_begin, HBTInt &task_end);
extern void logspace(double xmin,double xmax,int N, vector <float> &x);

#ifdef HAS_GSL
extern void EigenAxis(double Ixx, double Ixy, double Ixz, double Iyy, double Iyz, double Izz, float Axis[3][3]);
#endif

#endif
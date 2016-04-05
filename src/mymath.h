#ifndef MYMATH_HEADER_INCLUDED
#define MYMATH_HEADER_INCLUDED

#include <iostream>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <sys/stat.h>
#include <chrono>
#include <array>
#include <vector>
#include <iterator>

#include "datatypes.h"
#include "config_parser.h"

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
  void Reset()
  {
	tickers.clear();
  }
  size_t Size()
  {
	return tickers.size();
  }
  double GetSeconds(int i_start)
  {
	if(i_start+1<Size())
	  return chrono::duration_cast<chrono::duration<double> >(tickers[i_start+1]-tickers[i_start]).count();
	else
	  return 0;
  }
  double GetSeconds(int i_start, int i_end)
  {
	return chrono::duration_cast<chrono::duration<double> >(tickers[i_end]-tickers[i_start]).count();
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
extern std::ostream& operator << (std::ostream& o, HBTxyz &a);
extern void swap_Nbyte(void *data2swap,size_t nel,size_t mbyte);
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
inline void copyHBTxyz(HBTxyz & dest, const HBTxyz & src)
{
  memcpy(dest, src, sizeof(HBTxyz));
}

template <class T, std::size_t N>
ostream& operator<<(ostream& o, const array<T, N>& arr)
{
  o<<"(";
  copy(arr.cbegin(), arr.cend(), ostream_iterator<T>(o, ", "));
  o<<")";
  return o;
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
	HBTReal dx[3];
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];
	dx[2]=x[2]-y[2];
	return dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
}
inline HBTReal Distance(const HBTReal x[3], const HBTReal y[3])
{
	return sqrt(Distance2(x,y));
}
#define NEAREST(x) (((x)>HBTConfig.BoxHalf)?((x)-HBTConfig.BoxSize):(((x)<-HBTConfig.BoxHalf)?((x)+HBTConfig.BoxSize):(x)))
inline HBTReal PeriodicDistance(const HBTReal x[3], const HBTReal y[3])
{
	HBTReal dx[3];
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];
	dx[2]=x[2]-y[2];
	if(HBTConfig.PeriodicBoundaryOn)
	{
	  dx[0]=NEAREST(dx[0]);
	  dx[1]=NEAREST(dx[1]);
	  dx[2]=NEAREST(dx[2]);
	}
	return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
}

#ifdef HAS_GSL
extern void EigenAxis(double Ixx, double Ixy, double Ixz, double Iyy, double Iyz, double Izz, float Axis[3][3]);
#endif

#endif
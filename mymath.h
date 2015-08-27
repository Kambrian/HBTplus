#ifndef MYMATH_HEADER_INCLUDED
#define MYMATH_HEADER_INCLUDED

#include <iostream>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <sys/stat.h>

#include "datatypes.h"

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
inline HBTReal position_modulus(HBTReal x, HBTReal boxsize)
{//shift the positions to within [0,boxsize)
	HBTReal y;
	if(x>=0&&x<boxsize) return x;
	y=x/boxsize;
	return (y-floor(y))*boxsize;
}
inline HBTReal distance(HBTReal x[3],HBTReal y[3])
{
	HBTReal dx[3];
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];
	dx[2]=x[2]-y[2];
	#ifdef PERIODIC_BDR
	dx[0]=NEAREST(dx[0]);
	dx[1]=NEAREST(dx[1]);
	dx[2]=NEAREST(dx[2]);
	#endif
	return sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
}
#endif
#ifndef MYMATH_HEADER_INCLUDED
#define MYMATH_HEADER_INCLUDED

#include <iostream>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "datatypes.h"

// #define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){ fprintf(logfile,"Error opening file '%s'\n",filename);	fflush(logfile); exit(1);}
// #ifdef PERIODIC_BDR
// #define NEAREST(x) (((x)>BOXHALF)?((x)-BOXSIZE):(((x)<-BOXHALF)?((x)+BOXSIZE):(x)))
			/*this macro can well manipulate boundary condition because 
			* usually a halo is much smaller than boxhalf
			* so that any distance within the halo should be smaller than boxhalf */
// #endif
#define get_bit(x,k) (((x)&(1<<k))>>k)

extern std::ostream& operator << (std::ostream& o, HBTxyz &a);
extern void swap_Nbyte(void *data2swap,size_t nel,size_t mbyte);
extern size_t fread_swap(void *buf,size_t Nsize,size_t Nbuf,FILE *fp, int FlagByteSwap);

#endif
#ifndef FORTREAD_HEADER_INCLUDED
#define FORTREAD_HEADER_INCLUDED

extern "C" 
{
extern int alloc_file_unit_(int *reset);
extern void open_fortran_file_(const char *filename,int *fileno,int *endian,int *error);
extern void skip_fortran_record_(int *fileno);
extern void read_fortran_record1_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record2_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record4_(void *arr,long int *arr_len,int *fileno);
extern void read_fortran_record8_(void *arr,long int *arr_len,int *fileno);
extern void close_fortran_file_(int *fileno);
extern void read_part_arr_imajor_(float partarr[][3],long int *np,int *fileno);//old format
extern void read_part_arr_xmajor_(float partarr[][3],long int *np,int *fileno);//newest format
extern void read_part_header_int4_(int *np,int *ips,float *ztp,float *omgt,float *lbdt,
			      float *boxsize,float *xscale,float *vscale,int *fileno);
extern void read_part_header_int8_(long int *np,long int *ips,float *ztp,float *omgt,float *lbdt,
			      float *boxsize,float *xscale,float *vscale,int *fileno);
extern void read_ic_header_int4_(int *np,int *ips,float *ztp,float *omgt,float *lbdt,
			      float *boxsize,float *xscale,float *vscale,int *fileno);
extern void read_ic_header_int8_(long int *np,long int *ips,float *ztp,float *omgt,float *lbdt,
			      float *boxsize,float *xscale,float *vscale,int *fileno);
extern void read_group_header_int4_(float *b,int *ngrp, int *fileno);
extern void read_group_header_int8_(float *b,long int *ngrp, int *fileno);
}
//Caution:: only support SAME_INTTYPE and SAME_REALTYPE; if they differ, need further modification of code!
#ifdef HBT_INT8
#define read_fortran_record_HBTInt read_fortran_record8_
#define read_part_header read_part_header_int8_
#define read_ic_header read_ic_header_int8_
#define read_group_header read_group_header_int8_
#else
#define read_fortran_record_HBTInt read_fortran_record4_
#define read_part_header read_part_header_int4_
#define read_ic_header read_ic_header_int4_
#define read_group_header read_group_header_int4_
#endif

#endif
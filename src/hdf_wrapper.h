#ifndef HDF_WRAPPER_INCLUDED
#define HDF_WRAPPER_INCLUDED

#include "hdf5.h"
// #include "hdf5_hl.h"	
// #include "H5Cpp.h"
#include <iostream>

#ifdef HBT_REAL8
#define H5T_HBTReal H5T_NATIVE_DOUBLE
#else
#define H5T_HBTReal H5T_NATIVE_FLOAT
#endif
#ifdef HBT_INT8
#define H5T_HBTInt H5T_NATIVE_LONG
#else 
#define H5T_HBTInt H5T_NATIVE_INT
#endif

extern void writeHDFmatrix(hid_t file, const void * buf, const char * name, hsize_t ndim, const hsize_t *dims, hid_t dtype, hid_t dtype_file);
extern herr_t SetStringAttribute(hid_t loc_id, const char *obj_name, const char *attr_name, const char *content);

inline int GetDatasetDims(hid_t dset, hsize_t dims[])
{
  hid_t dspace=H5Dget_space(dset);
  int ndim=H5Sget_simple_extent_dims(dspace, dims, NULL);
  H5Sclose(dspace);
  return ndim;
}
inline herr_t ReclaimVlenData(hid_t dset, hid_t dtype, void * buf)
{
  herr_t status;
  hid_t dspace=H5Dget_space(dset);
  status=H5Dvlen_reclaim(dtype, dspace, H5P_DEFAULT, buf);
  status=H5Sclose(dspace);
  return status;
}
inline herr_t ReadDataset(hid_t file, const char *name, hid_t dtype, void *buf)
/* read named dataset from file into buf.
 * dtype specifies the datatype of buf; it does not need to be the same as the storage type in file*/
{
  herr_t status;
  hid_t dset=H5Dopen2(file, name, H5P_DEFAULT);
  status=H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  if(status<0)
  {
    const int bufsize=1024;
    char grpname[bufsize],filename[bufsize];
    H5Iget_name(file, grpname, bufsize);
    H5Fget_name(file, filename, bufsize);
    std::cerr<<"####ERROR READING "<<grpname<<"/"<<name<<" from "<<filename<<", error number "<<status<<std::endl<<std::flush;
  }
  H5Dclose(dset);
  return status;
}
inline herr_t ReadAttribute(hid_t loc_id, const char *obj_name, const char *attr_name, hid_t dtype, void *buf)
/* read named attribute of object into buf. if loc_id fully specifies the object, obj_name="."
 * dtype specifies the datatype of buf; it does not need to be the same as the storage type in file*/
{
  herr_t status;
  hid_t attr=H5Aopen_by_name(loc_id, obj_name, attr_name, H5P_DEFAULT, H5P_DEFAULT);
  status=H5Aread(attr, dtype, buf);
  status=H5Aclose(attr);
  return status;
}
inline void writeHDFmatrix(hid_t file, const void * buf, const char * name, hsize_t ndim, const hsize_t *dims, hid_t dtype)
{
  writeHDFmatrix(file, buf, name, ndim, dims, dtype, dtype);
}
#endif

#include "hdf_wrapper.h"

void writeHDFmatrix(hid_t file, const void * buf, const char * name, hsize_t ndim, const hsize_t *dims, hid_t dtype, hid_t dtype_file)
{
  hid_t dataspace = H5Screate_simple (ndim, dims, NULL);
  hid_t dataset= H5Dcreate2(file, name, dtype_file, dataspace, H5P_DEFAULT, H5P_DEFAULT,
                H5P_DEFAULT);
  if(!(NULL==buf||0==dims[0]))
  {
	herr_t status = H5Dwrite (dataset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	if(status<0)
	{
	  const int bufsize=1024;
	  char grpname[bufsize],filename[bufsize];
	  H5Iget_name(file, grpname, bufsize);
	  H5Fget_name(file, filename, bufsize);
	  std::cerr<<"####ERROR WRITING "<<grpname<<"/"<<name<<" into "<<filename<<", error number "<<status<<std::endl<<std::flush;
	}
  }
  H5Dclose(dataset);
  H5Sclose(dataspace);
}
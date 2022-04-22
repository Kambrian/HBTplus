#include "hdf_wrapper.h"
#include <cstring>

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

herr_t SetStringAttribute(hid_t loc_id, const char *obj_name, const char *attr_name, const char *content)
/* set string attribute named attr_name to object obj_name at loc_id. 
 * if loc_id fully specifies the object, obj_name="."
 * content specifies the attribute content
 *
 * equivalent to H5LTset_attribute_string() function in H5LT
 */
{
     hid_t dataspace  = H5Screate(H5S_SCALAR);
     hid_t attr_type = H5Tcopy(H5T_C_S1);
     H5Tset_size(attr_type, strlen(content)+1);
     H5Tset_strpad(attr_type, H5T_STR_NULLTERM);
     hid_t attr = H5Acreate_by_name(loc_id, obj_name, attr_name, attr_type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     herr_t status = H5Awrite(attr, attr_type, content);
     
     H5Sclose(dataspace);
     H5Tclose(attr_type);
     H5Aclose(attr);
  
  return status;
}

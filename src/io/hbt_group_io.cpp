/* IO for groups created by HBT's own FoF code.
 * The reader will be called in halo_io.cpp under "my_group_format".
 * 
 * To use this reader, set
 * 
 *   GroupFileFormat hbt_group_format
 * 
 * in config file.
 */

#include "../mymath.h"
#include "hbt_group_io.h"
#include "../hdf_wrapper.h"

namespace HBTGroupIO
{
  void GetHDFFileName(string &filename, int SnapshotIndex)
  {
    stringstream formater;
    formater<<HBTConfig.SubhaloPath<<"/"<<"HaloSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
    filename=formater.str();
  }
  
  HBTInt LoadHDFGroups(int snapshot_index, int snapshot_id, vector <Halo_t> &Halos)
  {
    string filename;
    GetHDFFileName(filename, snapshot_index);
    hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    {//check snapid
      int snapid;
      ReadDataset(file, "SnapshotId", H5T_NATIVE_INT, &snapid);
      assert(snapid==snapshot_id);
    }
    
    
    hid_t dset=H5Dopen2(file, "HaloParticles", H5P_DEFAULT);
    hsize_t dims[1];
    GetDatasetDims(dset, dims);
    HBTInt nhalos=dims[0];
    vector <hvl_t> vl(nhalos);
    HBTInt ntotal=0;
    if(nhalos)
    {
      hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
      H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
      H5Tclose(H5T_HBTIntArr);
      for(HBTInt i=0;i<nhalos;i++)
      {
	Halos[i].Particles.resize(vl[i].len);
	memcpy(Halos[i].Particles.data(), vl[i].p, sizeof(HBTInt)*vl[i].len);
	ntotal+=vl[i].len;
      }
      ReclaimVlenData(dset, H5T_HBTIntArr, vl.data());
      H5Tclose(H5T_HBTIntArr);
    }
    H5Dclose(dset);
    
    H5Fclose(file);
    
    return ntotal;
  }
  
  void SaveHDFGroups(int snapshot_index, int snapshot_id, vector <vector<HBTInt> > &halos)
  {
    string filename;
    GetHDFFileName(filename, snapshot_index);
    cout<<"Saving "<<halos.size()<<" halos to "<<filename<<"..."<<endl;
    hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    hsize_t ndim=1, dim_atom[]={1};
    writeHDFmatrix(file, &snapshot_id, "SnapshotId", ndim, dim_atom, H5T_NATIVE_INT);
    HBTInt Ngroups=halos.size();
    hsize_t dim_grp[]={Ngroups};
    
    #ifdef UNSIGNED_LONG_ID_OUTPUT  
    hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_NATIVE_ULONG);//this does not affect anything inside the code, but the presentation in the hdf file
    #else
    hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
    #endif
    vector <hvl_t> vl(Ngroups);
    for(HBTInt i=0;i<Ngroups;i++)
    {
      vl[i].len=halos[i].size();
      vl[i].p=halos[i].data();
    }
    writeHDFmatrix(file, vl.data(), "HaloParticles", ndim, dim_grp, H5T_HBTIntArr);
    H5LTset_attribute_string(file,"HaloParticles","Comment","List of particle ids in each group");
    
    H5Fclose(file);
    H5Tclose(H5T_HBTIntArr);
  }
}
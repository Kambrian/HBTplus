#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "../datatypes.h"
#include "../snapshot_number.h"
#include "../subhalo.h"

void SubhaloSnapshot_t::BuildHDFDataType()
{
  H5T_SubhaloInMem=H5Tcreate(H5T_COMPOUND, sizeof (Subhalo_t));
  hsize_t dims=3;
  hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, &dims);
  #define InsertMember(x,t) H5Tinsert(H5T_SubhaloInMem, #x, HOFFSET(Subhalo_t, x), t)//;cout<<#x<<": "<<HOFFSET(Subhalo_t, x)<<endl
  InsertMember(TrackId, H5T_HBTInt);
  InsertMember(Nbound, H5T_HBTInt);
  InsertMember(HostHaloId, H5T_HBTInt);
  InsertMember(Rank, H5T_HBTInt);
  InsertMember(ComovingPosition, H5T_HBTxyz);
  InsertMember(PhysicalVelocity, H5T_HBTxyz);
  InsertMember(LastMaxMass, H5T_HBTInt);
  InsertMember(SnapshotIndexOfLastMaxMass, H5T_HBTInt);
  InsertMember(SnapshotIndexOfLastIsolation, H5T_HBTInt);
  #undef InsertMember	
  H5T_SubhaloInDisk=H5Tcopy(H5T_SubhaloInMem);
  H5Tpack(H5T_SubhaloInDisk); //clear fields not added.
//   Subhalo_t s;
//   cout<<(char *)&s.TrackId-(char *)&s<<","<<(char *)&s.Nbound-(char *)&s<<","<<(char *)&s.ComovingPosition-(char *)&s<<","<<(char *)&s.Particles-(char *)&s<<endl;
  /*  
  #define InsertMember(x,t) H5T_ParticleInMem.insertMember(#x, HOFFSET(Subhalo_t, x), t)//;cout<<#x<<": "<<HOFFSET(Subhalo_t, x)<<endl
  InsertMember(Id, H5T_HBTInt);
//   InsertMember(Mass, H5T_HBTReal);
//   InsertMember(ComovingPosition, H5T_HBTxyz);
//   InsertMember(PhysicalVelocity, H5T_HBTxyz);
  #undef InsertMember	
  H5T_ParticleInDisk.copy(H5T_ParticleInMem);
  H5T_ParticleInDisk.pack(); //clear fields not added.  
*/
  H5Tclose(H5T_HBTxyz);
}
void SubhaloSnapshot_t::GetSubFileName(string &filename)
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SubSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
  filename=formater.str();
}
void SubhaloSnapshot_t::GetSrcFileName(string &filename)
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SrcSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
  filename=formater.str();
}
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
  status=H5Dclose(dset);
  return status;
}
void SubhaloSnapshot_t::Load(int snapshot_index, bool load_src)
{
  if(snapshot_index<HBTConfig.MinSnapshotIndex)
  {
	cout<<"Warning: snapshot index "<<snapshot_index<<" is below MinSnapshotIndex of "<<HBTConfig.MinSnapshotIndex;
	cout<<", I wil not load anything. (Ignore this message if you are starting HBT from MinSnapshotIndex.)\n";
	return;
  }
  
  SetSnapshotIndex(snapshot_index);
  string filename;
  GetSubFileName(filename);
  hid_t dset, file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HBTInt snapshot_id;
  ReadDataset(file, "SnapshotId", H5T_HBTInt, &snapshot_id);
  assert(snapshot_id==SnapshotId);
  
  hsize_t dims[1];
  dset=H5Dopen2(file, "Subhalos", H5P_DEFAULT);
  GetDatasetDims(dset, dims);
  HBTInt nsubhalos=dims[0];
  Subhalos.resize(nsubhalos);
  if(nsubhalos)	H5Dread(dset, H5T_SubhaloInMem, H5S_ALL, H5S_ALL, H5P_DEFAULT, Subhalos.data());
  H5Dclose(dset);
 
  ReadDataset(file, "/Membership/NumberOfNewSubhalos", H5T_HBTInt, &MemberTable.NBirth);
  ReadDataset(file, "/Membership/NumberOfFakeHalos", H5T_HBTInt, &MemberTable.NFake);
  dset=H5Dopen2(file, "/Membership/GroupedTrackIds", H5P_DEFAULT);
  GetDatasetDims(dset, dims);
  HBTInt nhalos=dims[0]-1;
  vector <hvl_t> vl(dims[0]);
  hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
  if(nsubhalos)
  {
	H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
	MemberTable.Init(nhalos, nsubhalos);
	#define FILL_SUBGROUP(vl_id, halo_id, offset) {\
	MemberTable.SubGroups[halo_id].Bind(vl[vl_id].len, &(MemberTable.AllMembers[offset])); \
	memcpy(MemberTable.SubGroups[halo_id].data(), vl[vl_id].p, sizeof(HBTInt)*vl[vl_id].len);\
	offset+=vl[vl_id].len;}
	HBTInt offset=0;
	FILL_SUBGROUP(nhalos, -1, offset);
	for(HBTInt i=0;i<nhalos;i++)
	  FILL_SUBGROUP(i, i, offset);
	#undef FILL_SUBGROUP
	ReclaimVlenData(dset, H5T_HBTIntArr, vl.data());
	assert(offset==nsubhalos);
  }
  H5Dclose(dset);
  
  if(0==nsubhalos) return;

  vl.resize(nsubhalos);
  if(!load_src)
  {
	dset=H5Dopen2(file, "SubhaloParticles", H5P_DEFAULT);
	GetDatasetDims(dset, dims);
	assert(dims[0]==nsubhalos);
	H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
	for(HBTInt i=0;i<nsubhalos;i++)
	{
	  Subhalos[i].Particles.resize(vl[i].len);
	  memcpy(Subhalos[i].Particles.data(), vl[i].p, sizeof(HBTInt)*vl[i].len);
	}
	ReclaimVlenData(dset, H5T_HBTIntArr, vl.data());
	H5Dclose(dset);
  }
  else
  {
	GetSrcFileName(filename);
	hid_t file2=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	dset=H5Dopen2(file2,"SrchaloParticles", H5P_DEFAULT);
	GetDatasetDims(dset, dims);
	assert(dims[0]==nsubhalos);
	H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
	for(HBTInt i=0;i<nsubhalos;i++)
	{
	  Subhalos[i].Particles.resize(vl[i].len);
	  memcpy(Subhalos[i].Particles.data(), vl[i].p, sizeof(HBTInt)*vl[i].len);
	}
	ReclaimVlenData(dset, H5T_HBTIntArr, vl.data());
	H5Dclose(dset);
	H5Fclose(file2);
  }
  H5Fclose(file);
  H5Tclose(H5T_HBTIntArr);
  cout<<Subhalos.size()<<" subhaloes loaded at snapshot "<<SnapshotIndex<<"("<<SnapshotId<<")\n";
}

void writeHDFmatrix(hid_t file, const void * buf, const char * name, hsize_t ndim, const hsize_t *dims, hid_t dtype, hid_t dtype_file)
{
  hid_t dataspace = H5Screate_simple (ndim, dims, NULL);
  hid_t dataset= H5Dcreate2(file, name, dtype_file, dataspace, H5P_DEFAULT, H5P_DEFAULT,
                H5P_DEFAULT);
  if(!(NULL==buf||0==dims[0]))
  {
	herr_t status = H5Dwrite (dataset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  }
  H5Sclose(dataspace);
  H5Dclose(dataset);
}
inline void writeHDFmatrix(hid_t file, const void * buf, const char * name, hsize_t ndim, const hsize_t *dims, hid_t dtype)
{
  writeHDFmatrix(file, buf, name, ndim, dims, dtype, dtype);
}
void SubhaloSnapshot_t::Save()
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SubSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
  string filename(formater.str());
  cout<<"Saving to "<<filename<<"..."<<endl;
  hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t ndim=1, dim_atom[]={1};
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5T_NATIVE_INT);
  writeHDFmatrix(file, &Hz, "HubbleParam", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(file, &ScaleFactor, "ScaleFactor", ndim, dim_atom, H5T_HBTReal);
//   writeHDFmatrix(file, &MemberTable.NBirth, "NumberOfNewSubhalos", ndim, dim_atom, H5T_HBTInt);
//   writeHDFmatrix(file, &MemberTable.NFake, "NumberOfFakeHalos", ndim, dim_atom, H5T_HBTInt);
  
  hsize_t dim_sub[]={Subhalos.size()};
  writeHDFmatrix(file, Subhalos.data(), "Subhalos", ndim, dim_sub, H5T_SubhaloInMem, H5T_SubhaloInDisk); 
  
  hid_t datagrp=H5Gcreate2(file, "/Membership", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_string(file,"/Membership","Comment","List of subhaloes in each group.");
  writeHDFmatrix(datagrp, &MemberTable.NBirth, "NumberOfNewSubhalos", ndim, dim_atom, H5T_HBTInt);
  writeHDFmatrix(datagrp, &MemberTable.NFake, "NumberOfFakeHalos", ndim, dim_atom, H5T_HBTInt);
 
  MemberTable.SubIdToTrackId(Subhalos);
  HBTInt Ngroups=MemberTable.SubGroups.size();
  hsize_t dim_grp[]={(hsize_t)Ngroups+1};
  hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
  vector <hvl_t> vl(Ngroups+1);
  for(HBTInt i=0;i<Ngroups;i++)
  {
	vl[i].len=MemberTable.SubGroups[i].size();
	vl[i].p=MemberTable.SubGroups[i].data();
  }
  vl[Ngroups].len=MemberTable.SubGroups[-1].size();
  vl[Ngroups].p=MemberTable.SubGroups[-1].data();
  writeHDFmatrix(datagrp, vl.data(), "GroupedTrackIds", ndim, dim_grp, H5T_HBTIntArr);
  H5LTset_attribute_string(datagrp.getId(),"GroupedTrackIds","Comment","Nhalo+1 groups. The last group contain tracks outside any host halo (i.e., field subhaloes).");
  H5Gclose(datagrp);
  
  //now write the particle list for each subhalo
  vl.resize(Subhalos.size());
  for(HBTInt i=0;i<vl.size();i++)
  {
	vl[i].len=Subhalos[i].Nbound;
	vl[i].p=Subhalos[i].Particles.data();
  }
  writeHDFmatrix(file, vl.data(), "SubhaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  
  H5Fclose(file);
  
  formater.str("");
  formater.clear();
  formater<<HBTConfig.SubhaloPath<<"/SrcSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
  filename=formater.str();
  cout<<"Saving to "<<filename<<"..."<<endl;
  file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5T_NATIVE_INT);
  for(HBTInt i=0;i<vl.size();i++)
	vl[i].len=Subhalos[i].Particles.size();
  writeHDFmatrix(file, vl.data(), "SrchaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  H5Fclose(file);
  H5Tclose(H5T_HBTIntArr);
  cout<<Subhalos.size()<<" subhaloes saved: "<<MemberTable.NBirth<<" birth, "<< MemberTable.NFake<<" fake.\n";
}


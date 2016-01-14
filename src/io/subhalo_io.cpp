#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "../mpi_wrapper.h"
#include "../datatypes.h"
#include "../snapshot_number.h"
#include "../subhalo.h"


void SubhaloSnapshot_t::BuildHDFDataType()
{
  H5T_SubhaloInMem=H5Tcreate(H5T_COMPOUND, sizeof (Subhalo_t));
  hsize_t dims[2]={3,3};
  hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);
  hid_t H5T_FloatVec3=H5Tarray_create2(H5T_NATIVE_FLOAT, 1, dims);
  #define InsertMember(x,t) H5Tinsert(H5T_SubhaloInMem, #x, HOFFSET(Subhalo_t, x), t)//;cout<<#x<<": "<<HOFFSET(Subhalo_t, x)<<endl
  InsertMember(TrackId, H5T_HBTInt);
  InsertMember(Nbound, H5T_HBTInt);
  InsertMember(HostHaloId, H5T_HBTInt);
  InsertMember(Rank, H5T_HBTInt);
  InsertMember(LastMaxMass, H5T_HBTInt);  
  InsertMember(SnapshotIndexOfLastMaxMass, H5T_NATIVE_INT);
  InsertMember(SnapshotIndexOfLastIsolation, H5T_NATIVE_INT);
  InsertMember(SnapshotIndexOfBirth, H5T_NATIVE_INT);
  InsertMember(SnapshotIndexOfDeath, H5T_NATIVE_INT);
  InsertMember(RmaxComoving, H5T_NATIVE_FLOAT);
  InsertMember(VmaxPhysical, H5T_NATIVE_FLOAT);
  InsertMember(LastMaxVmaxPhysical, H5T_NATIVE_FLOAT);
  InsertMember(SnapshotIndexOfLastMaxVmax, H5T_NATIVE_INT);
  InsertMember(R2SigmaComoving, H5T_NATIVE_FLOAT);
  InsertMember(RHalfComoving, H5T_NATIVE_FLOAT);
  InsertMember(R200CritComoving, H5T_NATIVE_FLOAT);
  InsertMember(R200MeanComoving, H5T_NATIVE_FLOAT);
  InsertMember(RVirComoving, H5T_NATIVE_FLOAT);
  InsertMember(M200Crit, H5T_NATIVE_FLOAT);
  InsertMember(M200Mean, H5T_NATIVE_FLOAT);
  InsertMember(MVir, H5T_NATIVE_FLOAT);
  InsertMember(SpecificSelfPotentialEnergy, H5T_NATIVE_FLOAT);
  InsertMember(SpecificSelfKineticEnergy, H5T_NATIVE_FLOAT);
  InsertMember(SpecificAngularMomentum, H5T_FloatVec3);
#ifdef ENABLE_EXPERIMENTAL_PROPERTIES
  InsertMember(SpinPeebles, H5T_FloatVec3);
  InsertMember(SpinBullock, H5T_FloatVec3);
#endif
#ifdef HAS_GSL
  dims[0]=3;
  dims[1]=3;
  hid_t H5T_FloatVec33=H5Tarray_create2(H5T_NATIVE_FLOAT, 2, dims);
  InsertMember(InertialEigenVector, H5T_FloatVec33);
  InsertMember(InertialEigenVectorWeighted, H5T_FloatVec33);
  H5Tclose(H5T_FloatVec33);
#endif
  dims[0]=6;
  hid_t H5T_FloatVec6=H5Tarray_create2(H5T_NATIVE_FLOAT, 1, dims);
  InsertMember(InertialTensor,H5T_FloatVec6);
  InsertMember(InertialTensorWeighted, H5T_FloatVec6);
  H5Tclose(H5T_FloatVec6);

  InsertMember(ComovingMostBoundPosition, H5T_HBTxyz);
  InsertMember(PhysicalMostBoundVelocity, H5T_HBTxyz);
  InsertMember(ComovingAveragePosition, H5T_HBTxyz);
  InsertMember(PhysicalAverageVelocity, H5T_HBTxyz);
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
  H5Tclose(H5T_FloatVec3);
  H5Tclose(H5T_HBTxyz);
}
void SubhaloSnapshot_t::GetSubFileName(string &filename, int iFile)
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SubSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<"."<<iFile<<".hdf5"; //or use snapshotid
  filename=formater.str();
}
void SubhaloSnapshot_t::GetSrcFileName(string &filename, int iFile)
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SrcSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<"."<<iFile<<".hdf5"; //or use snapshotid
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
void SubhaloSnapshot_t::Load(MpiWorker_t &world, int snapshot_index, bool load_src)
{
  if(snapshot_index<HBTConfig.MinSnapshotIndex)
  {
	if(world.rank()==0)
	  cout<<"Skipping empty snapshot "<<snapshot_index<<"\n";
	return;
  }
  
  int iFile=world.rank();
  
  SetSnapshotIndex(snapshot_index);
  string filename;
  GetSubFileName(filename, iFile);
  hid_t dset, file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HBTInt snapshot_id;
  ReadDataset(file, "SnapshotId", H5T_HBTInt, &snapshot_id);
  assert(snapshot_id==SnapshotId);
  
  ReadDataset(file, "OmegaM0", H5T_HBTReal, &OmegaM0);
  ReadDataset(file, "OmegaLambda0", H5T_HBTReal, &OmegaLambda0);
  ReadDataset(file, "HubbleParam", H5T_HBTReal, &Hz);
  ReadDataset(file, "ScaleFactor", H5T_HBTReal, &ScaleFactor);
  
  int NumberOfFiles;
  ReadDataset(file, "NumberOfFiles", H5T_HBTInt, &NumberOfFiles);
  if(NumberOfFiles!=world.size())
  {
	cerr<<"Error: can only read subhalo files with the same number of processes as saved\n";
	exit(1);
  }
  
//   ReadDataset(file, "NumberOfNewSubhalos", H5T_HBTInt, &MemberTable.NBirth);
//   ReadDataset(file, "NumberOfFakeHalos", H5T_HBTInt, &MemberTable.NFake);
  
  hsize_t dims[1];
  dset=H5Dopen2(file, "Subhalos", H5P_DEFAULT);
  GetDatasetDims(dset, dims);
  HBTInt nsubhalos=dims[0];
  Subhalos.resize(nsubhalos);
  if(nsubhalos)	H5Dread(dset, H5T_SubhaloInMem, H5S_ALL, H5S_ALL, H5P_DEFAULT, Subhalos.data());
  H5Dclose(dset);
 
  if(0==nsubhalos) return;
  
  vector <hvl_t> vl(dims[0]);
  vl.resize(nsubhalos);
  hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
  if(!load_src)
  {
	dset=H5Dopen2(file, "SubhaloParticles", H5P_DEFAULT);
	GetDatasetDims(dset, dims);
	assert(dims[0]==nsubhalos);
	H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
	for(HBTInt i=0;i<nsubhalos;i++)
	{
	  Subhalos[i].Particles.resize(vl[i].len);
	  HBTInt *p=(HBTInt *)(vl[i].p);
	  for(HBTInt j=0;j<vl[i].len;j++)
		Subhalos[i].Particles[j].Id=p[j];
	}
	ReclaimVlenData(dset, H5T_HBTIntArr, vl.data());
	H5Dclose(dset);
  }
  else
  {
	GetSrcFileName(filename, iFile);
	hid_t file2=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	dset=H5Dopen2(file2,"SrchaloParticles", H5P_DEFAULT);
	GetDatasetDims(dset, dims);
	assert(dims[0]==nsubhalos);
	H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
	for(HBTInt i=0;i<nsubhalos;i++)
	{
	  Subhalos[i].Particles.resize(vl[i].len);
	  HBTInt *p=(HBTInt *)(vl[i].p);
	  for(HBTInt j=0;j<vl[i].len;j++)
		Subhalos[i].Particles[j].Id=p[j];
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
void SubhaloSnapshot_t::Save(MpiWorker_t &world)
{
  int iFile=world.rank();
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SubSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<"."<<iFile<<".hdf5"; //or use snapshotid
  string filename(formater.str());
  cout<<"Saving "<<Subhalos.size()<<" subhaloes to "<<filename<<"..."<<endl;
  hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t ndim=1, dim_atom[]={1};
  int nfiles=world.size();
  writeHDFmatrix(file, &nfiles, "NumberOfFiles", ndim, dim_atom, H5T_NATIVE_INT);
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5T_NATIVE_INT);
  writeHDFmatrix(file, &OmegaM0, "OmegaM0", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(file, &OmegaLambda0, "OmegaLambda0", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(file, &Hz, "HubbleParam", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(file, &ScaleFactor, "ScaleFactor", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(file, &MemberTable.NBirth, "NumberOfNewSubhalos", ndim, dim_atom, H5T_HBTInt);
  writeHDFmatrix(file, &MemberTable.NFake, "NumberOfFakeHalos", ndim, dim_atom, H5T_HBTInt);
  
  hsize_t dim_sub[]={Subhalos.size()};
  writeHDFmatrix(file, Subhalos.data(), "Subhalos", ndim, dim_sub, H5T_SubhaloInMem, H5T_SubhaloInDisk); 
  
  /* do not write membertable; this can always be rebuilt.
  H5::Group datagrp(file.createGroup("/Membership"));
  H5LTset_attribute_string(file.getId(),"/Membership","Comment","List of subhaloes in each group.");
  writeHDFmatrix(datagrp, &MemberTable.NBirth, "NumberOfNewSubhalos", ndim, dim_atom, H5T_HBTInt);
  writeHDFmatrix(datagrp, &MemberTable.NFake, "NumberOfFakeHalos", ndim, dim_atom, H5T_HBTInt);
 
//   MemberTable.SubIdToTrackId(Subhalos);
  HBTInt Ngroups=MemberTable.SubGroups.size();
  hsize_t dim_grp[]={(hsize_t)Ngroups+1};
  H5::VarLenType H5T_HBTIntArr(&H5T_HBTInt);
  vector <hvl_t> vl(Ngroups+1);
  for(HBTInt i=0;i<Ngroups;i++)
  {
	vl[i].len=MemberTable.SubGroups[i].size();
	vl[i].p=MemberTable.SubGroups[i].data();
  }
  vl[Ngroups].len=MemberTable.SubGroups[-1].size();
  vl[Ngroups].p=MemberTable.SubGroups[-1].data();
  writeHDFmatrix(datagrp, vl.data(), "GroupedSubIds", ndim, dim_grp, H5T_HBTIntArr);
  H5LTset_attribute_string(datagrp.getId(),"GroupedSubIds","Comment","Nhalo+1 groups. The last group contains tracks outside any host halo (i.e., field subhaloes). Each group contains the list of subhaloes inside this group");
  */
  //now write the particle list for each subhalo
  hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
  vector <hvl_t> vl(Subhalos.size());
  vector <HBTInt> IdBuffer;
  {
	HBTInt NumberOfParticles=0;
	for(HBTInt i=0;i<Subhalos.size();i++)
	  NumberOfParticles+=Subhalos[i].Particles.size();
	IdBuffer.reserve(NumberOfParticles);
	HBTInt offset=0;
	for(HBTInt i=0;i<Subhalos.size();i++)
	{
	  vl[i].len=Subhalos[i].Nbound;
	  vl[i].p=&IdBuffer[offset];
	  offset+=Subhalos[i].Particles.size();
	  for(auto && p: Subhalos[i].Particles)
		IdBuffer.push_back(p.Id);
	}
  }
  writeHDFmatrix(file, vl.data(), "SubhaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  
  H5Fclose(file);
  
  formater.str("");
  formater.clear();
  formater<<HBTConfig.SubhaloPath<<"/SrcSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<"."<<iFile<<".hdf5"; //or use snapshotid
  filename=formater.str();
  cout<<"Saving "<<Subhalos.size()<<" subhaloes to "<<filename<<"..."<<endl;
  file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5T_NATIVE_INT);
  for(HBTInt i=0;i<vl.size();i++)
	vl[i].len=Subhalos[i].Particles.size();
  writeHDFmatrix(file, vl.data(), "SrchaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  H5Fclose(file);
  H5Tclose(H5T_HBTIntArr);
  cout<<Subhalos.size()<<" subhaloes saved: "<<MemberTable.NBirth<<" birth, "<< MemberTable.NFake<<" fake.\n";
}


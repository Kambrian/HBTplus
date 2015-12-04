#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "../datatypes.h"
#include "../snapshot_number.h"
#include "../subhalo.h"
#include "../boost_mpi.h"

void SubhaloSnapshot_t::BuildHDFDataType()
{
  hsize_t dims=3;
  H5::ArrayType H5T_HBTxyz(H5T_HBTReal, 1, &dims);
  #define InsertMember(x,t) H5T_SubhaloInMem.insertMember(#x, HOFFSET(Subhalo_t, x), t)//;cout<<#x<<": "<<HOFFSET(Subhalo_t, x)<<endl
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
  H5T_SubhaloInDisk.copy(H5T_SubhaloInMem);
  H5T_SubhaloInDisk.pack(); //clear fields not added.
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
void SubhaloSnapshot_t::Load(mpi::communicator &world, int snapshot_index, bool load_src)
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
  H5::H5File file(filename.c_str(), H5F_ACC_RDONLY);
  
  H5::DataSet dset(file.openDataSet("SnapshotId"));
  HBTInt snapshot_id;
  dset.read(&snapshot_id, H5T_HBTInt);//HBTInt does not have to be the same as the datatype in the file.
  dset.close();
  assert(snapshot_id==SnapshotId);
  
  dset=file.openDataSet("NumberOfFiles");
  int NumberOfFiles;
  dset.read(&NumberOfFiles, H5::PredType::NATIVE_INT);
  dset.close();
  if(NumberOfFiles!=world.size())
  {
	cerr<<"Error: can only read subhalo files with the same number of processes as saved\n";
	exit(1);
  }
  
  hsize_t dims[1];
  dset=file.openDataSet("Subhalos");
  dset.getSpace().getSimpleExtentDims(dims);
  HBTInt nsubhalos=dims[0];
  Subhalos.resize(nsubhalos);
  if(nsubhalos)	dset.read(Subhalos.data(), H5T_SubhaloInMem);
  dset.close();
 
  
  if(0==nsubhalos) return;
  
  vector <hvl_t> vl(dims[0]);
  H5::VarLenType H5T_HBTIntArr(&H5T_HBTInt);
  vl.resize(nsubhalos);
  if(!load_src)
  {
	dset=file.openDataSet("SubhaloParticles");
	dset.getSpace().getSimpleExtentDims(dims);
	assert(dims[0]==nsubhalos);
	dset.read(vl.data(), H5T_HBTIntArr);
	for(HBTInt i=0;i<nsubhalos;i++)
	{
	  Subhalos[i].Particles.resize(vl[i].len);
	  HBTInt *p=(HBTInt *)(vl[i].p);
	  for(HBTInt j=0;j<vl[i].len;j++)
		Subhalos[i].Particles[j].Id=p[j];
	}
	dset.vlenReclaim(vl.data(), H5T_HBTIntArr,dset.getSpace());
	dset.close();
  }
  else
  {
	GetSrcFileName(filename, iFile);
	H5::H5File file2(filename.c_str(), H5F_ACC_RDONLY);
	dset=file2.openDataSet("SrchaloParticles");
	dset.getSpace().getSimpleExtentDims(dims);
	assert(dims[0]==nsubhalos);
	dset.read(vl.data(), H5T_HBTIntArr);
	for(HBTInt i=0;i<nsubhalos;i++)
	{
	  Subhalos[i].Particles.resize(vl[i].len);
	  HBTInt *p=(HBTInt *)(vl[i].p);
	  for(HBTInt j=0;j<vl[i].len;j++)
		Subhalos[i].Particles[j].Id=p[j];
	}
	dset.vlenReclaim(vl.data(), H5T_HBTIntArr,dset.getSpace());
	dset.close();
  }
  cout<<Subhalos.size()<<" subhaloes loaded at snapshot "<<SnapshotIndex<<"("<<SnapshotId<<")\n";
}

void writeHDFmatrix(H5::CommonFG &file, const void * buf, const char * name, const hsize_t ndim, const hsize_t *dims, const H5::DataType &dtype, const H5::DataType &dtype_file)
{
  H5::DataSpace dataspace(ndim, dims);
  H5::DataSet dataset(file.createDataSet( name, dtype_file, dataspace));
  if(NULL==buf||0==dims[0]) return;
  dataset.write(buf, dtype);
}
inline void writeHDFmatrix(H5::CommonFG &file, const void * buf, const char * name, const hsize_t ndim, const hsize_t *dims, const H5::DataType &dtype)
{
  writeHDFmatrix(file, buf, name, ndim, dims, dtype, dtype);
}
void SubhaloSnapshot_t::Save(mpi::communicator &world)
{
  int iFile=world.rank();
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SubSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<"."<<iFile<<".hdf5"; //or use snapshotid
  string filename(formater.str());
  cout<<"Saving "<<Subhalos.size()<<" subhaloes to "<<filename<<"..."<<endl;
  H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);

  hsize_t ndim=1, dim_atom[]={1};
  int nfiles=world.size();
  writeHDFmatrix(file, &nfiles, "NumberOfFiles", ndim, dim_atom, H5::PredType::NATIVE_INT);
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5::PredType::NATIVE_INT);
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
  H5::VarLenType H5T_HBTIntArr(&H5T_HBTInt);
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
  
  file.close();
  
  formater.str("");
  formater.clear();
  formater<<HBTConfig.SubhaloPath<<"/SrcSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<"."<<iFile<<".hdf5"; //or use snapshotid
  filename=formater.str();
  cout<<"Saving "<<Subhalos.size()<<" subhaloes to "<<filename<<"..."<<endl;
  file=H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5::PredType::NATIVE_INT);
  for(HBTInt i=0;i<vl.size();i++)
	vl[i].len=Subhalos[i].Particles.size();
  writeHDFmatrix(file, vl.data(), "SrchaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  file.close();
  cout<<Subhalos.size()<<" subhaloes saved: "<<MemberTable.NBirth<<" birth, "<< MemberTable.NFake<<" fake.\n";
}


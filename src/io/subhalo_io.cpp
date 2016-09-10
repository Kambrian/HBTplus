#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "../datatypes.h"
#include "../snapshot_number.h"
#include "../subhalo.h"
#include "../hdf_wrapper.h"

SubhaloSnapshot_t::~SubhaloSnapshot_t()
{
  H5Tclose(H5T_SubhaloInDisk);
  H5Tclose(H5T_SubhaloInMem);
}

void SubhaloSnapshot_t::BuildHDFDataType()
{
  H5T_SubhaloInMem=H5Tcreate(H5T_COMPOUND, sizeof (Subhalo_t));
  hsize_t dims[2]={3,3};
  hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);
  hid_t H5T_FloatVec3=H5Tarray_create2(H5T_NATIVE_FLOAT, 1, dims);
  #define InsertMember(x,t) H5Tinsert(H5T_SubhaloInMem, #x, HOFFSET(Subhalo_t, x), t)//;cout<<#x<<": "<<HOFFSET(Subhalo_t, x)<<endl
  InsertMember(TrackId, H5T_HBTInt);
  InsertMember(Nbound, H5T_HBTInt);
  InsertMember(Mbound, H5T_NATIVE_FLOAT);
  
  dims[0]=TypeMax;
  hid_t H5T_HBTIntArray_TypeMax=H5Tarray_create2(H5T_HBTInt, 1, dims);
  hid_t H5T_FloatArray_TypeMax=H5Tarray_create2(H5T_NATIVE_FLOAT, 1, dims);
#ifndef DM_ONLY
  InsertMember(NboundType, H5T_HBTIntArray_TypeMax);
  InsertMember(MboundType, H5T_FloatArray_TypeMax);
#endif
  H5Tclose(H5T_HBTIntArray_TypeMax);
  H5Tclose(H5T_FloatArray_TypeMax);
  
  InsertMember(HostHaloId, H5T_HBTInt);
  InsertMember(Rank, H5T_HBTInt);
  InsertMember(LastMaxMass, H5T_NATIVE_FLOAT);  
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
  InsertMember(SpinPeebles, H5T_FloatVec3);
  InsertMember(SpinBullock, H5T_FloatVec3);
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
  
  ReadDataset(file, "/Cosmology/OmegaM0", H5T_HBTReal, &Cosmology.OmegaM0);
  ReadDataset(file, "/Cosmology/OmegaLambda0", H5T_HBTReal, &Cosmology.OmegaLambda0);
  ReadDataset(file, "/Cosmology/HubbleParam", H5T_HBTReal, &Cosmology.Hz);
  ReadDataset(file, "/Cosmology/ScaleFactor", H5T_HBTReal, &Cosmology.ScaleFactor);
#ifdef DM_ONLY
  ReadDataset(file, "/Cosmology/ParticleMass", H5T_HBTReal, &Cosmology.ParticleMass);
#endif
  
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
  
  {//read nested subhalos
    dset=H5Dopen2(file, "NestedSubhalos", H5P_DEFAULT);
    GetDatasetDims(dset, dims);
    assert(dims[0]==nsubhalos);
    H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
    for(HBTInt i=0;i<nsubhalos;i++)
    {
      Subhalos[i].NestedSubhalos.resize(vl[i].len);
      memcpy(Subhalos[i].NestedSubhalos.data(), vl[i].p, sizeof(HBTInt)*vl[i].len);
    }
    ReclaimVlenData(dset, H5T_HBTIntArr, vl.data());
    H5Dclose(dset);
  }
  
#ifdef SAVE_BINDING_ENERGY
  {//read binding energies
	dset=H5Dopen2(file, "BindingEnergies", H5P_DEFAULT);
    GetDatasetDims(dset, dims);
    assert(dims[0]==nsubhalos);
	hid_t H5T_FloatArr=H5Tvlen_create(H5T_NATIVE_FLOAT);
    H5Dread(dset, H5T_FloatArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
    for(HBTInt i=0;i<nsubhalos;i++)
    {
      Subhalos[i].Energies.resize(vl[i].len);
      memcpy(Subhalos[i].Energies.data(), vl[i].p, sizeof(float)*vl[i].len);
    }
    ReclaimVlenData(dset, H5T_FloatArr, vl.data());
	H5Tclose(H5T_FloatArr);
    H5Dclose(dset);
  }
#endif
  
  H5Fclose(file);
  H5Tclose(H5T_HBTIntArr);
  cout<<Subhalos.size()<<" subhaloes loaded at snapshot "<<SnapshotIndex<<"("<<SnapshotId<<")\n";
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
  hid_t cosmology=H5Gcreate2(file, "/Cosmology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  writeHDFmatrix(cosmology, &Cosmology.OmegaM0, "OmegaM0", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(cosmology, &Cosmology.OmegaLambda0, "OmegaLambda0", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(cosmology, &Cosmology.Hz, "HubbleParam", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(cosmology, &Cosmology.ScaleFactor, "ScaleFactor", ndim, dim_atom, H5T_HBTReal);
#ifdef DM_ONLY
  writeHDFmatrix(cosmology, &Cosmology.ParticleMass, "ParticleMass", ndim, dim_atom, H5T_HBTReal);
#endif
  H5Gclose(cosmology);
//   writeHDFmatrix(file, &MemberTable.NBirth, "NumberOfNewSubhalos", ndim, dim_atom, H5T_HBTInt);
//   writeHDFmatrix(file, &MemberTable.NFake, "NumberOfFakeHalos", ndim, dim_atom, H5T_HBTInt);
  
  hid_t datagrp=H5Gcreate2(file, "/Membership", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_string(file,"/Membership","Comment","List of subhaloes in each group.");
  writeHDFmatrix(datagrp, &MemberTable.NBirth, "NumberOfNewSubhalos", ndim, dim_atom, H5T_HBTInt);
  writeHDFmatrix(datagrp, &MemberTable.NFake, "NumberOfFakeHalos", ndim, dim_atom, H5T_HBTInt);
 
  MemberTable.SubIdToTrackId(Subhalos);
  HBTInt Ngroups=MemberTable.SubGroups.size();
  hsize_t dim_grp[]={(hsize_t)Ngroups+1};

#ifdef UNSIGNED_LONG_ID_OUTPUT  
  hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_NATIVE_ULONG);//this does not affect anything inside the code, but the presentation in the hdf file
#else
  hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
#endif
  vector <hvl_t> vl(Ngroups+1);
  for(HBTInt i=0;i<Ngroups;i++)
  {
	vl[i].len=MemberTable.SubGroups[i].size();
	vl[i].p=MemberTable.SubGroups[i].data();
  }
  vl[Ngroups].len=MemberTable.SubGroups[-1].size();
  vl[Ngroups].p=MemberTable.SubGroups[-1].data();
  writeHDFmatrix(datagrp, vl.data(), "GroupedTrackIds", ndim, dim_grp, H5T_HBTIntArr);
  H5LTset_attribute_string(datagrp,"GroupedTrackIds","Comment","Nhalo+1 groups. The last group contain tracks outside any host halo (i.e., field subhaloes).");
  H5Gclose(datagrp);
   
  hsize_t dim_sub[]={Subhalos.size()};
  //now write the particle list for each subhalo
  if(HBTConfig.SaveSubParticleProperties)
  {
	struct Particle_t
	{
	  HBTInt ParticleIndex;
	  HBTxyz ComovingPosition;
	  HBTxyz PhysicalVelocity;
#ifndef DM_ONLY
	  HBTReal Mass;
#ifdef HAS_THERMAL_ENERGY
	  HBTReal InternalEnergy;
#endif
	  int Type;
#endif
	};
	hid_t H5T_Particle=H5Tcreate(H5T_COMPOUND, sizeof (Particle_t));
	hsize_t dim_xyz=3;
	hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, &dim_xyz);
	#define InsertMember(x,t) H5Tinsert(H5T_Particle, #x, HOFFSET(Particle_t, x), t)
	InsertMember(ParticleIndex, H5T_HBTInt);
	InsertMember(ComovingPosition, H5T_HBTxyz);
	InsertMember(PhysicalVelocity, H5T_HBTxyz);
#ifndef DM_ONLY
	InsertMember(Mass, H5T_HBTReal);
	#ifdef HAS_THERMAL_ENERGY
	InsertMember(InternalEnergy, H5T_HBTReal);
	#endif
	InsertMember(Type, H5T_NATIVE_INT);
#endif
	#undef InsertMember
	H5Tclose(H5T_HBTxyz);
	
	hid_t H5T_ParticleArr=H5Tvlen_create(H5T_Particle);

	vector <vector<Particle_t> > SubParticles(Subhalos.size());
	vl.resize(Subhalos.size());
	for(HBTInt subid=0;subid<Subhalos.size();subid++)
	{
	  auto &IdList=Subhalos[subid].Particles;
	  HBTInt np=Subhalos[subid].Nbound;
	  auto &ParticleList=SubParticles[subid];
	  ParticleList.resize(np);
	  for(HBTInt i=0;i<np;i++)
	  {
		auto &p=ParticleList[i];
		auto &ind=IdList[i];
		p.ParticleIndex=ind;
		copyHBTxyz(p.ComovingPosition, SnapshotPointer->GetComovingPosition(ind));
		copyHBTxyz(p.PhysicalVelocity, SnapshotPointer->GetPhysicalVelocity(ind));
#ifndef DM_ONLY
		p.Mass=SnapshotPointer->GetMass(ind);
		#ifdef HAS_THERMAL_ENERGY
		p.InternalEnergy=SnapshotPointer->GetInternalEnergy(ind);
		#endif
		p.Type=SnapshotPointer->GetParticleType(ind);
#endif
	  }
	  vl[subid].len=np;
	  vl[subid].p=ParticleList.data();
	}
	writeHDFmatrix(file, vl.data(), "ParticleProperties", ndim, dim_sub, H5T_ParticleArr);
	H5Tclose(H5T_ParticleArr);
	H5Tclose(H5T_Particle);
  }
  
  vl.resize(Subhalos.size());
  for(HBTInt i=0;i<vl.size();i++)
  {
	vl[i].len=Subhalos[i].NestedSubhalos.size();
	vl[i].p=Subhalos[i].NestedSubhalos.data();
  }
  writeHDFmatrix(file, vl.data(), "NestedSubhalos", ndim, dim_sub, H5T_HBTIntArr);
  H5LTset_attribute_string(file,"/NestedSubhalos","Comment","List of the indices of first-level sub-subhaloes within each subhalo.");
  
  ParticleIndexToId();
  for(HBTInt i=0;i<vl.size();i++)
  {
	vl[i].len=Subhalos[i].Nbound;
	vl[i].p=Subhalos[i].Particles.data();
  }
  writeHDFmatrix(file, vl.data(), "SubhaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  
#ifdef SAVE_BINDING_ENERGY
  hid_t H5T_FloatArr=H5Tvlen_create(H5T_NATIVE_FLOAT);
  for(HBTInt i=0;i<vl.size();i++)
	vl[i].p=Subhalos[i].Energies.data();
  writeHDFmatrix(file, vl.data(), "BindingEnergies", ndim, dim_sub, H5T_FloatArr);
  H5Tclose(H5T_FloatArr);
#endif
  
  writeHDFmatrix(file, Subhalos.data(), "Subhalos", ndim, dim_sub, H5T_SubhaloInMem, H5T_SubhaloInDisk); //write after calling ParticleIndexToId, for MostBoundParticleId
  
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


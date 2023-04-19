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
  InsertMember(Depth, H5T_NATIVE_INT);
  InsertMember(LastMaxMass, H5T_NATIVE_FLOAT);
  InsertMember(SnapshotIndexOfLastMaxMass, H5T_NATIVE_INT);
  InsertMember(SnapshotIndexOfLastIsolation, H5T_NATIVE_INT);
  InsertMember(SnapshotIndexOfBirth, H5T_NATIVE_INT);
  InsertMember(SnapshotIndexOfDeath, H5T_NATIVE_INT);
  InsertMember(SnapshotIndexOfSink, H5T_NATIVE_INT);
  InsertMember(RmaxComoving, H5T_NATIVE_FLOAT);
  InsertMember(VmaxPhysical, H5T_NATIVE_FLOAT);
  InsertMember(LastMaxVmaxPhysical, H5T_NATIVE_FLOAT);
  InsertMember(SnapshotIndexOfLastMaxVmax, H5T_NATIVE_INT);
  InsertMember(R2SigmaComoving, H5T_NATIVE_FLOAT);
  InsertMember(RHalfComoving, H5T_NATIVE_FLOAT);
  InsertMember(BoundR200CritComoving, H5T_NATIVE_FLOAT);
//   InsertMember(R200MeanComoving, H5T_NATIVE_FLOAT);
//   InsertMember(RVirComoving, H5T_NATIVE_FLOAT);
  InsertMember(BoundM200Crit, H5T_NATIVE_FLOAT);
//   InsertMember(M200Mean, H5T_NATIVE_FLOAT);
//   InsertMember(MVir, H5T_NATIVE_FLOAT);
  InsertMember(SpecificSelfPotentialEnergy, H5T_NATIVE_FLOAT);
  InsertMember(SpecificSelfKineticEnergy, H5T_NATIVE_FLOAT);
  InsertMember(SpecificAngularMomentum, H5T_FloatVec3);
//   InsertMember(SpinPeebles, H5T_FloatVec3);
//   InsertMember(SpinBullock, H5T_FloatVec3);
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

  InsertMember(ComovingAveragePosition, H5T_HBTxyz);
  InsertMember(PhysicalAverageVelocity, H5T_HBTxyz);
  InsertMember(ComovingMostBoundPosition, H5T_HBTxyz);
  InsertMember(PhysicalMostBoundVelocity, H5T_HBTxyz);
  InsertMember(MostBoundParticleId, H5T_HBTInt);

  InsertMember(SinkTrackId, H5T_HBTInt);
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
inline void Subhalo_t::DuplicateMostBoundParticleId()
{
  if(Particles.size())
    MostBoundParticleId=Particles[0].Id;
  else
    MostBoundParticleId=SpecialConst::NullParticleId;
}
string SubhaloSnapshot_t::GetSubDir()
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/"<<setw(3)<<setfill('0')<<SnapshotIndex;
  return formater.str();
}
void SubhaloSnapshot_t::GetSubFileName(string &filename, int iFile, const string &ftype)
{
  stringstream formater;
  formater<<GetSubDir()<<"/"+ftype+"Snap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<"."<<iFile<<".hdf5"; //or use snapshotid
  filename=formater.str();
}
void SubhaloSnapshot_t::Load(MpiWorker_t &world, int snapshot_index, const SubReaderDepth_t depth)
{
  if(snapshot_index<HBTConfig.MinSnapshotIndex)
  {
	if(world.rank()==0)
	  cout<<"Skipping empty snapshot "<<snapshot_index<<"\n";
	return;
  }
  SetSnapshotIndex(snapshot_index);

  int NumberOfFiles;
  HBTInt TotNumberOfSubs;
  if(world.rank()==0)
  {
  string filename;
  GetSubFileName(filename, 0);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  ReadDataset(file, "NumberOfFiles", H5T_NATIVE_INT, &NumberOfFiles);
  ReadDataset(file, "NumberOfSubhalosInAllFiles", H5T_HBTInt, &TotNumberOfSubs);
  }
  MPI_Bcast(&NumberOfFiles, 1, MPI_INT, 0, world.Communicator);

  Subhalos.clear();
  HBTInt nfiles_skip, nfiles_end;
  AssignTasks(world.rank(), world.size(), NumberOfFiles, nfiles_skip, nfiles_end);

  for(int i=0, ireader=0;i<world.size();i++, ireader++)
  {
	if(ireader==HBTConfig.MaxConcurrentIO)
	{
	  ireader=0;//reset reader count
	  MPI_Barrier(world.Communicator);//wait for every thread to arrive.
	}
	if(i==world.rank())//read
	{
	  for(int iFile=nfiles_skip;iFile<nfiles_end;iFile++)
		ReadFile(iFile, depth);
	}
  }

//   cout<<Subhalos.size()<<" subhaloes loaded at snapshot "<<SnapshotIndex<<"("<<SnapshotId<<")\n";

  HBTInt NumSubs=Subhalos.size(), NumSubsAll_loaded=0;
  MPI_Reduce(&NumSubs, &NumSubsAll_loaded, 1, MPI_HBT_INT, MPI_SUM, 0, world.Communicator);
  if(world.rank()==0)
  {
	if(NumSubsAll_loaded!=TotNumberOfSubs)
	{
	  ostringstream msg;
	  msg<<"Error reading SubSnap "<<snapshot_index<<": total number of subhaloes expected="<<TotNumberOfSubs<<", loaded="<<NumSubsAll_loaded<<endl;
	  throw runtime_error(msg.str().c_str());
	}
	cout<<TotNumberOfSubs<<" subhalos loaded at snapshot "<<SnapshotIndex<<"("<<SnapshotId<<")\n";
  }
}
void SubhaloSnapshot_t::ReadFile(int iFile, const SubReaderDepth_t depth)
{//Read iFile for current snapshot.

  string filename;
  GetSubFileName(filename, iFile);
  hid_t dset, file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HBTInt snapshot_id;
  ReadDataset(file, "SnapshotId", H5T_HBTInt, &snapshot_id);
  assert(snapshot_id==SnapshotId);

  ReadDataset(file, "/Cosmology/OmegaM0", H5T_HBTReal, &Cosmology.OmegaM0);
  ReadDataset(file, "/Cosmology/OmegaLambda0", H5T_HBTReal, &Cosmology.OmegaLambda0);
  ReadDataset(file, "/Cosmology/HubbleParam", H5T_HBTReal, &Cosmology.Hz);
  ReadDataset(file, "/Cosmology/ScaleFactor", H5T_HBTReal, &Cosmology.ScaleFactor);
  Cosmology.Set(Cosmology.ScaleFactor, Cosmology.OmegaM0, Cosmology.OmegaLambda0);

//   ReadDataset(file, "NumberOfNewSubhalos", H5T_HBTInt, &MemberTable.NBirth);
//   ReadDataset(file, "NumberOfFakeHalos", H5T_HBTInt, &MemberTable.NFake);

  hsize_t dims[1];
  dset=H5Dopen2(file, "Subhalos", H5P_DEFAULT);
  GetDatasetDims(dset, dims);
  HBTInt nsubhalos=dims[0], nsubhalos_old=Subhalos.size();
  Subhalos.resize(nsubhalos+nsubhalos_old);
  if(nsubhalos)	H5Dread(dset, H5T_SubhaloInMem, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Subhalos[nsubhalos_old]);
  H5Dclose(dset);

  if(nsubhalos)
  {
    Subhalo_t * NewSubhalos=&Subhalos[nsubhalos_old];
    vector <hvl_t> vl(dims[0]);
    vl.resize(nsubhalos);
    hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
    if(depth==SubReaderDepth_t::SubParticles||depth==SubReaderDepth_t::SrcParticles)
    {
      hid_t file2;
      switch(depth)
      {
	case SubReaderDepth_t::SubParticles:
	  dset=H5Dopen2(file, "SubhaloParticles", H5P_DEFAULT);
	  break;
	case SubReaderDepth_t::SrcParticles:
	  GetSubFileName(filename, iFile, "Src");
	  file2=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	  dset=H5Dopen2(file2,"SrchaloParticles", H5P_DEFAULT);
	  break;
      }
      GetDatasetDims(dset, dims);
      assert(dims[0]==nsubhalos);
      H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
      for(HBTInt i=0;i<nsubhalos;i++)
      {
	NewSubhalos[i].Particles.resize(vl[i].len);
	HBTInt *p=(HBTInt *)(vl[i].p);
	for(HBTInt j=0;j<vl[i].len;j++)
	  NewSubhalos[i].Particles[j].Id=p[j];
      }
      ReclaimVlenData(dset, H5T_HBTIntArr, vl.data());
      H5Dclose(dset);
      if(depth==SubReaderDepth_t::SrcParticles)
	H5Fclose(file2);
    }

    {//read nested subhalos
      dset=H5Dopen2(file, "NestedSubhalos", H5P_DEFAULT);
      GetDatasetDims(dset, dims);
      assert(dims[0]==nsubhalos);
      H5Dread(dset, H5T_HBTIntArr, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl.data());
      for(HBTInt i=0;i<nsubhalos;i++)
      {
	NewSubhalos[i].NestedSubhalos.resize(vl[i].len);
	memcpy(NewSubhalos[i].NestedSubhalos.data(), vl[i].p, sizeof(HBTInt)*vl[i].len);
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
	NewSubhalos[i].Energies.resize(vl[i].len);
	memcpy(NewSubhalos[i].Energies.data(), vl[i].p, sizeof(float)*vl[i].len);
      }
      ReclaimVlenData(dset, H5T_FloatArr, vl.data());
	  H5Tclose(H5T_FloatArr);
      H5Dclose(dset);
    }
  #endif

    H5Tclose(H5T_HBTIntArr);
  }

  H5Fclose(file);
}

void SubhaloSnapshot_t::Save(MpiWorker_t &world)
{
  string subdir=GetSubDir();
  mkdir(subdir.c_str(), 0755);

  HBTInt NumSubsAll=0, NumSubs=Subhalos.size();
  MPI_Allreduce(&NumSubs, &NumSubsAll, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);

  HBTInt MaxSize=0, MaxSizeAll=0;
  for(HBTInt i=0;i<Subhalos.size();i++)
  {
    if(MaxSize<Subhalos[i].Nbound)
      MaxSize=Subhalos[i].Nbound;
  }
  MPI_Reduce(&MaxSize, &MaxSizeAll, 1, MPI_HBT_INT, MPI_MAX, 0, world.Communicator);
  if(world.rank()==0) cout<<"Max subhalo size="<<MaxSizeAll<<endl;

  if(world.rank()==0) cout<<"saving "<<NumSubsAll<<" subhalos to "<<subdir<<endl;
  for(int i=0, ireader=0;i<world.size();i++, ireader++)
  {
	if(ireader==HBTConfig.MaxConcurrentIO)
	{
	  ireader=0;//reset reader count
	  MPI_Barrier(world.Communicator);//wait for every thread to arrive.
	}
	if(i==world.rank())//read
	{
	  WriteFile(world.rank(), world.size(), NumSubsAll);
	}
  }
}
void SubhaloSnapshot_t::WriteFile(int iFile, int nfiles, HBTInt NumSubsAll)
{
  string filename;
  GetSubFileName(filename, iFile);
//   cout<<"Saving "<<Subhalos.size()<<" subhaloes to "<<filename<<"..."<<endl;
  hid_t file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t ndim=1, dim_atom[]={1};
  writeHDFmatrix(file, &nfiles, "NumberOfFiles", ndim, dim_atom, H5T_NATIVE_INT);
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5T_NATIVE_INT);
  hid_t cosmology=H5Gcreate2(file, "/Cosmology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  writeHDFmatrix(cosmology, &Cosmology.OmegaM0, "OmegaM0", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(cosmology, &Cosmology.OmegaLambda0, "OmegaLambda0", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(cosmology, &Cosmology.Hz, "HubbleParam", ndim, dim_atom, H5T_HBTReal);
  writeHDFmatrix(cosmology, &Cosmology.ScaleFactor, "ScaleFactor", ndim, dim_atom, H5T_HBTReal);
  H5Gclose(cosmology);
  writeHDFmatrix(file, &MemberTable.NBirth, "NumberOfNewSubhalos", ndim, dim_atom, H5T_HBTInt);
  writeHDFmatrix(file, &MemberTable.NFake, "NumberOfFakeHalos", ndim, dim_atom, H5T_HBTInt);
  writeHDFmatrix(file, &NumSubsAll, "NumberOfSubhalosInAllFiles", ndim, dim_atom, H5T_HBTInt);//for data verification

  vector <hvl_t> vl(Subhalos.size());
  hsize_t dim_sub[]={Subhalos.size()};
    //now write the particle list for each subhalo
  if(HBTConfig.SaveSubParticleProperties)
  {
	hid_t H5T_ParticleInMem=H5Tcreate(H5T_COMPOUND, sizeof (Particle_t));
	hsize_t dim_xyz=3;
	hid_t H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, &dim_xyz);
	#define InsertMember(x,t) H5Tinsert(H5T_ParticleInMem, #x, HOFFSET(Particle_t, x), t)
	InsertMember(ComovingPosition, H5T_HBTxyz);
	InsertMember(PhysicalVelocity, H5T_HBTxyz);
	InsertMember(Mass, H5T_HBTReal);
#ifndef DM_ONLY
	#ifdef HAS_THERMAL_ENERGY
	InsertMember(InternalEnergy, H5T_HBTReal);
	#endif
	InsertMember(Type, H5T_NATIVE_INT);
#endif
	#undef InsertMember
	H5Tclose(H5T_HBTxyz);
	hid_t H5T_ParticleInDisk=H5Tcopy(H5T_ParticleInMem);
	H5Tpack(H5T_ParticleInDisk);
	hid_t H5T_ParticleArrInMem=H5Tvlen_create(H5T_ParticleInMem);
	hid_t H5T_ParticleArrInDisk=H5Tvlen_create(H5T_ParticleInDisk);

	for(HBTInt subid=0;subid<Subhalos.size();subid++)
	{
	  HBTInt np=Subhalos[subid].Nbound;
	  vl[subid].len=np;
	  vl[subid].p=Subhalos[subid].Particles.data();
	}
	writeHDFmatrix(file, vl.data(), "ParticleProperties", ndim, dim_sub, H5T_ParticleArrInMem, H5T_ParticleArrInDisk);
	H5Tclose(H5T_ParticleArrInMem);
	H5Tclose(H5T_ParticleArrInDisk);
	H5Tclose(H5T_ParticleInMem);
	H5Tclose(H5T_ParticleInDisk);
  }

#ifdef UNSIGNED_LONG_ID_OUTPUT
  hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_NATIVE_ULONG);//this does not affect anything inside the code, but the presentation in the hdf file
#else
  hid_t H5T_HBTIntArr=H5Tvlen_create(H5T_HBTInt);
#endif

  for(HBTInt i=0;i<vl.size();i++)
  {
	vl[i].len=Subhalos[i].NestedSubhalos.size();
	vl[i].p=Subhalos[i].NestedSubhalos.data();
  }
  writeHDFmatrix(file, vl.data(), "NestedSubhalos", ndim, dim_sub, H5T_HBTIntArr);
  H5LTset_attribute_string(file,"/NestedSubhalos","Comment","List of the TrackIds of first-level sub-subhaloes within each subhalo.");

  #ifdef SAVE_BINDING_ENERGY
  hid_t H5T_FloatArr=H5Tvlen_create(H5T_NATIVE_FLOAT);
  for(HBTInt i=0;i<vl.size();i++)
  {
	vl[i].len=Subhalos[i].Nbound;
	vl[i].p=Subhalos[i].Energies.data();
  }
  writeHDFmatrix(file, vl.data(), "BindingEnergies", ndim, dim_sub, H5T_FloatArr);
  H5Tclose(H5T_FloatArr);
  #endif

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
	  Subhalos[i].DuplicateMostBoundParticleId();//dump for galform
	}
  }
  writeHDFmatrix(file, vl.data(), "SubhaloParticles", ndim, dim_sub, H5T_HBTIntArr);

  writeHDFmatrix(file, Subhalos.data(), "Subhalos", ndim, dim_sub, H5T_SubhaloInMem, H5T_SubhaloInDisk);

  H5Fclose(file);

  GetSubFileName(filename, iFile, "Src");
//   cout<<"Saving "<<Subhalos.size()<<" subhaloes to "<<filename<<"..."<<endl;
  file=H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5T_NATIVE_INT);
  for(HBTInt i=0;i<vl.size();i++)
	vl[i].len=Subhalos[i].Particles.size();
  writeHDFmatrix(file, vl.data(), "SrchaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  H5Fclose(file);
  H5Tclose(H5T_HBTIntArr);
//   cout<<Subhalos.size()<<" subhaloes saved: "<<MemberTable.NBirth<<" birth, "<< MemberTable.NFake<<" fake.\n";
}


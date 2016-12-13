/*read subhalo snapshot from MPI output (multiple files per snapshot)
 * this file requires HDF5 library with multi-thread support; otherwise please comment out the #pragma omp... near line 76.
 */
#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "../datatypes.h"
#include "../snapshot_number.h"
#include "../subhalo.h"

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
HBTInt SubhaloSnapshot_t::GetNumberOfSubhalos(int iFile)
{
  string filename;
  GetSubFileName(filename, iFile);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);  
  hsize_t dims[1];
  hid_t dset=H5Dopen2(file, "Subhalos", H5P_DEFAULT);
  GetDatasetDims(dset, dims);
  H5Dclose(dset);
  H5Fclose(file);
  return dims[0];
}
void SubhaloSnapshot_t::LoadSubDir(int snapshot_index, const SubReaderDepth_t depth)
{
  if(snapshot_index<HBTConfig.MinSnapshotIndex)
	  cout<<"Skipping empty snapshot "<<snapshot_index<<"\n";
  SetSnapshotIndex(snapshot_index);
  
  int NumberOfFiles;
  HBTInt TotNumberOfSubs;

  string filename;
  GetSubFileName(filename, 0);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);  
  ReadDataset(file, "NumberOfFiles", H5T_NATIVE_INT, &NumberOfFiles);
  ReadDataset(file, "NumberOfSubhalosInAllFiles", H5T_HBTInt, &TotNumberOfSubs);
  H5Fclose(file);
  FileOffset.resize(NumberOfFiles);
  HBTInt offset=0;
  for(int iFile=0;iFile<NumberOfFiles;iFile++)
  {
    FileOffset[iFile]=offset;
    offset+=GetNumberOfSubhalos(iFile);
  }
  if(offset!=TotNumberOfSubs)
  {
    ostringstream msg;
    msg<<"Error reading SubSnap "<<snapshot_index<<": total number of subhaloes expected="<<TotNumberOfSubs<<", loaded="<<offset<<endl;
    throw runtime_error(msg.str().c_str());
  }
    
  Subhalos.clear();  
  Subhalos.resize(TotNumberOfSubs); 
#pragma omp parallel for
  for(int iFile=0;iFile<NumberOfFiles;iFile++)
      ReadFile(iFile, depth);
  
  cout<<Subhalos.size()<<" subhaloes loaded at snapshot "<<SnapshotIndex<<"("<<SnapshotId<<")\n";
  
  //now build membertable
  HBTInt nhost=0;
#pragma omp parallel for reduction(max:nhost)
  for(HBTInt i=0;i<TotNumberOfSubs;i++)
  {
    if(Subhalos[i].HostHaloId>nhost)
      nhost=Subhalos[i].HostHaloId;
  }
  nhost++;
  MemberTable.Build(nhost, Subhalos, true); 
}
void SubhaloSnapshot_t::ReadFile(int iFile, const SubReaderDepth_t depth)
{//Read iFile for current snapshot. 
  
  string filename;
  GetSubFileName(filename, iFile);
  hid_t dset, file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HBTInt snapshot_id;
  ReadDataset(file, "SnapshotId", H5T_HBTInt, &snapshot_id);
  assert(snapshot_id==SnapshotId);
  
  if(iFile==0)
  {
  ReadDataset(file, "/Cosmology/OmegaM0", H5T_HBTReal, &Cosmology.OmegaM0);
  ReadDataset(file, "/Cosmology/OmegaLambda0", H5T_HBTReal, &Cosmology.OmegaLambda0);
  ReadDataset(file, "/Cosmology/HubbleParam", H5T_HBTReal, &Cosmology.Hz);
  ReadDataset(file, "/Cosmology/ScaleFactor", H5T_HBTReal, &Cosmology.ScaleFactor);
  }
  
//   ReadDataset(file, "NumberOfNewSubhalos", H5T_HBTInt, &MemberTable.NBirth);
//   ReadDataset(file, "NumberOfFakeHalos", H5T_HBTInt, &MemberTable.NFake);
  
  hsize_t dims[1];
  dset=H5Dopen2(file, "Subhalos", H5P_DEFAULT);
  GetDatasetDims(dset, dims);
  HBTInt nsubhalos=dims[0];
  Subhalo_t * NewSubhalos=&Subhalos[FileOffset[iFile]];
  if(nsubhalos)	H5Dread(dset, H5T_SubhaloInMem, H5S_ALL, H5S_ALL, H5P_DEFAULT, NewSubhalos);
  H5Dclose(dset);
 
  if(0==nsubhalos) return;
  
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
	NewSubhalos[i].Particles[j]=p[j];
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
  
  H5Fclose(file);
  H5Tclose(H5T_HBTIntArr);
}


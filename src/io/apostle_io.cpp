using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <boost/concept_check.hpp>

#include "../snapshot.h"
#include "../mymath.h"
#include "../hdf_wrapper.h"
#include "apostle_io.h"

ApostleReader_t::ApostleReader_t(int snapshotId, vector <Particle_t> &particles, Cosmology_t &cosmology): SnapshotId(snapshotId), Particles(particles), Cosmology(cosmology)
{
  if(HBTConfig.SnapshotNameList.empty())
  {
	stringstream formatter;
	formatter<<"snapshot_"<<setw(3)<<setfill('0')<<SnapshotId;
	SnapshotName=formatter.str();
  }
  else
	SnapshotName=HBTConfig.SnapshotNameList[SnapshotId];
  
  hsize_t dims[]={TypeMax};
  MassTable_t=H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dims);
  CountTable_t=H5Tarray_create2(H5T_NATIVE_INT, 1, dims);
  dims[0]=3;
  H5T_HBTxyz=H5Tarray_create2(H5T_HBTReal, 1, dims);
}

ApostleReader_t::~ApostleReader_t()
{
  H5Tclose(MassTable_t);
  H5Tclose(CountTable_t);
  H5Tclose(H5T_HBTxyz);
}

void ApostleReader_t::GetFileName(int ifile, string &filename)
{
  string subname=SnapshotName;
  subname.erase(4, 4);//i.e., remove "shot" from "snapshot" or "snipshot"
  stringstream formatter;
  formatter<<HBTConfig.SnapshotPath<<"/"<<SnapshotName<<"/"<<subname<<"."<<ifile<<".hdf5";
  filename=formatter.str();
}

void ApostleReader_t::ReadHeader(int ifile, ApostleHeader_t &header)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  ReadAttribute(file, "Header", "NumFilesPerSnapshot", H5T_NATIVE_INT, &Header.NumberOfFiles);
  ReadAttribute(file, "Header", "BoxSize", H5T_NATIVE_DOUBLE, &Header.BoxSize);
  assert(Header.BoxSize==HBTConfig.BoxSize);
  ReadAttribute(file, "Header", "Time", H5T_NATIVE_DOUBLE, &Header.ScaleFactor);
  ReadAttribute(file, "Header", "Omega0", H5T_NATIVE_DOUBLE, &Header.OmegaM0);
  ReadAttribute(file, "Header", "OmegaLambda", H5T_NATIVE_DOUBLE, &Header.OmegaLambda0);  
  ReadAttribute(file, "Header", "MassTable", MassTable_t, Header.mass);
  ReadAttribute(file, "Header", "NumPart_ThisFile", CountTable_t, Header.npart);
  
  int np[TypeMax], np_high[TypeMax];
  ReadAttribute(file, "Header", "NumPart_Total", CountTable_t, np);
  ReadAttribute(file, "Header", "NumPart_Total_HighWord", CountTable_t, np_high);
  for(int i=0;i<TypeMax;i++)
	Header.npartTotal[i]=(((long)np_high[i])<<32)|np[i];
  
  H5Fclose(file);
}

HBTInt ApostleReader_t::CompileFileOffsets(int nfiles)
{
  HBTInt offset=0;
  np_file.reserve(nfiles);
  offset_file.reserve(nfiles);
  for(int ifile=0;ifile<nfiles;ifile++)
  {
	offset_file.push_back(offset);
	
	int np_this[TypeMax];
	string filename;
	GetFileName(ifile, filename);
	hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	ReadAttribute(file, "Header", "NumPart_ThisFile", CountTable_t, np_this);
	H5Fclose(file);
	HBTInt np=accumulate(begin(np_this), end(np_this), (HBTInt)0);
	
	np_file.push_back(np);
	offset+=np;
  }
  return offset;
}

void ApostleReader_t::ReadFile(int ifile)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  int np_this[TypeMax];
  ReadAttribute(file, "Header", "NumPart_ThisFile", CountTable_t, np_this);
 
  HBTReal vunit=sqrt(Header.ScaleFactor);
  HBTReal boxsize=Header.BoxSize;
  auto NewParticles=Particles.begin()+offset_file[ifile];
  for(int itype=0;itype<TypeMax;itype++)
  {
	int np=np_this[itype];
	stringstream grpname;
	grpname<<"PartType"<<itype;
	hid_t particle_data=H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);

	{//read position
	  vector <HBTxyz> x(np);
	  ReadDataset(particle_data, "Coordinates", H5T_HBTxyz, x.data());	
	  if(HBTConfig.PeriodicBoundaryOn)
	  {
		for(int i=0;i<np;i++)
		for(int j=0;j<3;j++)
		  x[i][j]=position_modulus(x[i][j], boxsize);
	  }
	  for(int i=0;i<np;i++)
		copyHBTxyz(NewParticles[i].ComovingPosition, x[i]);
	}
	
	{//velocity
	  vector <HBTxyz> v(np);
	  ReadDataset(particle_data, "Velocities", H5T_HBTxyz, v.data());
	  for(int i=0;i<np;i++)
		for(int j=0;j<3;j++)
		  NewParticles[i].PhysicalVelocity[j]=v[i][j]*vunit;
	}
	
	{//id
	  vector <HBTInt> id(np);
	  ReadDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data());
	  for(int i=0;i<np;i++)
		NewParticles[i].Id=id[i];
	}
	
	//mass
	if(Header.mass[itype]==0)
	{
	  vector <HBTReal> m(np);
	  ReadDataset(particle_data, "Masses", H5T_HBTReal, m.data());
	  for(int i=0;i<np;i++)
		NewParticles[i].Mass=m[i];
	}
	else
	{
	  for(int i=0;i<np;i++)
		NewParticles[i].Mass=Header.mass[itype];
	}
	
	//internal energy
	if(itype==0)
	{
	  vector <HBTReal> u(np);
	  ReadDataset(particle_data, "InternalEnergy", H5T_HBTReal, u.data());
	  for(int i=0;i<np;i++)
		NewParticles[i].InternalEnergy=u[i];
	}	
/*	else
	{
	  for(int i=0;i<np;i++)
		NewParticles[i].InternalEnergy=0; //necessary? maybe default initialized?
	}*/
	
	{//type
	  ParticleType_t t=static_cast<ParticleType_t>(itype);
	  for(int i=0;i<np;i++)
		NewParticles[i].Type=t;
	}
  }
  
  H5Fclose(file);
}

void ApostleReader_t::LoadSnapshot()
{
  ReadHeader(0, Header);
  Cosmology.Set(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);  
  HBTInt np=CompileFileOffsets(Header.NumberOfFiles);
  Particles.resize(np);
  
#pragma omp parallel for num_threads(HBTConfig.MaxConcurrentIO)
  for(int iFile=0; iFile<Header.NumberOfFiles; iFile++)
	ReadFile(iFile);
	
  cout<<" ( "<<Header.NumberOfFiles<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
}
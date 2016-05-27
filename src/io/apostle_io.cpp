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

void ApostleReader_t::SetSnapshot(int snapshotId)
{  
  if(HBTConfig.SnapshotNameList.empty())
  {
	stringstream formatter;
	formatter<<"snapshot_"<<setw(3)<<setfill('0')<<snapshotId;
	SnapshotName=formatter.str();
  }
  else
	SnapshotName=HBTConfig.SnapshotNameList[snapshotId];
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
  assert((HBTReal)Header.BoxSize==HBTConfig.BoxSize);
  ReadAttribute(file, "Header", "Time", H5T_NATIVE_DOUBLE, &Header.ScaleFactor);
  ReadAttribute(file, "Header", "Omega0", H5T_NATIVE_DOUBLE, &Header.OmegaM0);
  ReadAttribute(file, "Header", "OmegaLambda", H5T_NATIVE_DOUBLE, &Header.OmegaLambda0);  
  ReadAttribute(file, "Header", "MassTable", H5T_NATIVE_DOUBLE, Header.mass);
//   cout<<Header.mass[0]<<","<<Header.mass[1]<<","<<Header.mass[2]<<endl;
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, Header.npart);
  
  unsigned np[TypeMax], np_high[TypeMax];
  ReadAttribute(file, "Header", "NumPart_Total", H5T_NATIVE_UINT, np);
  ReadAttribute(file, "Header", "NumPart_Total_HighWord", H5T_NATIVE_UINT, np_high);
  for(int i=0;i<TypeMax;i++)
	Header.npartTotal[i]=(((unsigned long)np_high[i])<<32)|np[i];
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
	ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, np_this);
	H5Fclose(file);
	HBTInt np=accumulate(begin(np_this), end(np_this), (HBTInt)0);
	
	np_file.push_back(np);
	offset+=np;
  }
  return offset;
}

static void check_id_size(hid_t loc)
{
  hid_t dset=H5Dopen2(loc, "ParticleIDs", H5P_DEFAULT);
  hid_t dtype=H5Dget_type(dset);
  size_t ParticleIDStorageSize=H5Tget_size(dtype);
  assert(sizeof(HBTInt)>=ParticleIDStorageSize); //use HBTi8 or HBTdouble if you need long int for id
  H5Tclose(dtype);
  H5Dclose(dset);
}
void ApostleReader_t::ReadSnapshot(int ifile, Particle_t *ParticlesInFile)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  vector <int> np_this(TypeMax);
  vector <HBTInt> offset_this(TypeMax);
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, np_this.data());
  CompileOffsets(np_this, offset_this);
 
  HBTReal vunit=sqrt(Header.ScaleFactor);
  HBTReal boxsize=Header.BoxSize;
  for(int itype=0;itype<TypeMax;itype++)
  {
	int np=np_this[itype];
	auto ParticlesThisType=ParticlesInFile+offset_this[itype];
	stringstream grpname;
	grpname<<"PartType"<<itype;
	if(!H5Lexists(file, grpname.str().c_str(), H5P_DEFAULT)) continue;
	hid_t particle_data=H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);
// 	if(particle_data<0) continue; //skip non-existing type

	check_id_size(particle_data);
	
	{//read position
	  vector <HBTxyz> x(np);
	  ReadDataset(particle_data, "Coordinates", H5T_HBTReal, x.data());	
	  if(HBTConfig.PeriodicBoundaryOn)
	  {
		for(int i=0;i<np;i++)
		for(int j=0;j<3;j++)
		  x[i][j]=position_modulus(x[i][j], boxsize);
	  }
	  for(int i=0;i<np;i++)
		copyHBTxyz(ParticlesThisType[i].ComovingPosition, x[i]);
	}
	
	{//velocity
	  vector <HBTxyz> v(np);
	  ReadDataset(particle_data, "Velocities", H5T_HBTReal, v.data());
	  for(int i=0;i<np;i++)
		for(int j=0;j<3;j++)
		  ParticlesThisType[i].PhysicalVelocity[j]=v[i][j]*vunit;
	}
	
	{//id
	  vector <HBTInt> id(np);
	  ReadDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data());
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].Id=id[i];
	}
	
	//mass
	if(Header.mass[itype]==0)
	{
	  vector <HBTReal> m(np);
	  ReadDataset(particle_data, "Masses", H5T_HBTReal, m.data());
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].Mass=m[i];
	}
	else
	{
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].Mass=Header.mass[itype];
	}
	
	//internal energy
#ifdef UNBIND_WITH_THERMAL_ENERGY
	if(itype==0)
	{
	  vector <HBTReal> u(np);
	  ReadDataset(particle_data, "InternalEnergy", H5T_HBTReal, u.data());
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].InternalEnergy=u[i];
	}	
/*	else
	{
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].InternalEnergy=0; //necessary? maybe default initialized?
	}*/
#endif
	
	{//type
	  ParticleType_t t=static_cast<ParticleType_t>(itype);
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].Type=t;
	}
	
	H5Gclose(particle_data);
  }
  
  H5Fclose(file);
}

void ApostleReader_t::ReadGroupId(int ifile, ParticleHost_t *ParticlesInFile, bool FlagReadParticleId)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  vector <int> np_this(TypeMax);
  vector <HBTInt> offset_this(TypeMax);
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, np_this.data());
  CompileOffsets(np_this, offset_this);
  
  for(int itype=0;itype<TypeMax;itype++)
  {
	int np=np_this[itype];
	auto ParticlesThisType=ParticlesInFile+offset_this[itype];
	stringstream grpname;
	grpname<<"PartType"<<itype;
	if(!H5Lexists(file, grpname.str().c_str(), H5P_DEFAULT)) continue;
	hid_t particle_data=H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);
	
	if(FlagReadParticleId)
	{//id
	  check_id_size(particle_data);
	  vector <HBTInt> id(np);
	  ReadDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data());
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].ParticleId=id[i];
	}
	
	{//Hostid
	  vector <HBTInt> id(np);
	  ReadDataset(particle_data, "GroupNumber", H5T_HBTInt, id.data());
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].HostId=(id[i]<0?-id[i]:id[i]);//negative means unbound 
	}
	
	H5Gclose(particle_data);
  }
  
  H5Fclose(file);
}

void ApostleReader_t::LoadSnapshot(int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology)
{
  SetSnapshot(snapshotId);
  
  ReadHeader(0, Header);
  Cosmology.Set(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);  
  HBTInt np=CompileFileOffsets(Header.NumberOfFiles);
  Particles.resize(np);
  
#pragma omp parallel for num_threads(HBTConfig.MaxConcurrentIO)
  for(int iFile=0; iFile<Header.NumberOfFiles; iFile++)
  {
	ReadSnapshot(iFile, Particles.data()+offset_file[iFile]);
	cout<<iFile<<" "<<flush;
  }
  cout<<endl;
  cout<<" ( "<<Header.NumberOfFiles<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
//   cout<<" Particle[0]: x="<<Particles[0].ComovingPosition<<", v="<<Particles[0].PhysicalVelocity<<", m="<<Particles[0].Mass<<endl;
//   cout<<" Particle[2]: x="<<Particles[2].ComovingPosition<<", v="<<Particles[2].PhysicalVelocity<<", m="<<Particles[2].Mass<<endl;
//   cout<<" Particle[end]: x="<<Particles.back().ComovingPosition<<", v="<<Particles.back().PhysicalVelocity<<", m="<<Particles.back().Mass<<endl;
}

inline bool CompParticleHost(const ParticleHost_t &a, const ParticleHost_t &b)
{
  return a.HostId<b.HostId;
}

HBTInt ApostleReader_t::LoadGroups(int snapshotId, vector< Halo_t >& Halos)
{
  SetSnapshot(snapshotId);
  
  ReadHeader(0, Header);
  HBTInt NumberOfParticles=CompileFileOffsets(Header.NumberOfFiles);
  vector <ParticleHost_t> Particles(NumberOfParticles);
  bool FlagReadId=!HBTConfig.GroupLoadedIndex;
  
  cout<<"reading group files: ";
  #pragma omp parallel for num_threads(HBTConfig.MaxConcurrentIO)
  for(int iFile=0; iFile<Header.NumberOfFiles; iFile++)
  {
	ReadGroupId(iFile, Particles.data()+offset_file[iFile], FlagReadId);
	cout<<iFile<<" "<<flush;
  }
  cout<<endl;
  
  if(!FlagReadId)
  {//fill with index
	#pragma omp parallel for 
	for(HBTInt i=0;i<NumberOfParticles;i++)
	  Particles[i].ParticleId=i;
  }
  
  const int NullGroupId=1<<30; //1073741824
  sort(Particles.begin(), Particles.end(), CompParticleHost);
  assert(Particles.back().HostId==NullGroupId);//max haloid==NullGroupId
  assert(Particles.front().HostId>=0);//min haloid>=0
  
  HBTInt NumGroups=0;
  HBTInt imax=lower_bound(Particles.begin(), Particles.end(), Particles.back(), CompParticleHost)-Particles.begin();
  if(imax>0)
	NumGroups=Particles[imax-1].HostId+1;
  Halos.resize(NumGroups);
  
  vector <HBTInt> HaloOffset(NumGroups+1, 0);
  HaloOffset.back()=imax;//offset of NullGroup
  #pragma omp parallel
  {
	#pragma omp for
	for(HBTInt i=1;i<imax;i++)
	{
	  HBTInt hid1=Particles[i-1].HostId, hid2=Particles[i].HostId;
	  if(hid1!=hid2)
	  {
		for(HBTInt hid=hid1+1;hid<=hid2;hid++)
		  HaloOffset[hid]=i;
	  }
	}  
	#pragma omp for
	for(HBTInt i=0;i<NumGroups;i++)
	{
	  HBTInt np=HaloOffset[i+1]-HaloOffset[i];  
	  auto &p=Halos[i].Particles;
	  p.resize(np);
	  auto p_in=Particles.data()+HaloOffset[i];
	  for(HBTInt j=0;j<np;j++)
		p[j]=p_in[j].ParticleId;
	}
  }
  cout<<Halos.size()<<" groups loaded";
  if(Halos.size()) cout<<" : "<<Halos[0].Particles.size();
  if(Halos.size()>1) cout<<","<<Halos[1].Particles.size()<<"...";
  cout<<endl;
  
  return imax;
}

bool IsApostleGroup(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 7)=="apostle";
}
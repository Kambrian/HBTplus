using namespace std;
#include <numeric>
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <chrono>

#include "../snapshot.h"
#include "../mymath.h"
#include "../hdf_wrapper.h"
#include "apostle_io.h"

bool ApostleReader_t::IsIllustris(hid_t file)
{//TODO: rely on config param to decide
      stringstream grpname;
	grpname<<"PartType"<<3;
	if(H5Lexists(file, grpname.str().c_str(), H5P_DEFAULT)){
            hid_t type3_data=H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);
            if(H5Lexists(type3_data, "ParentID", H5P_DEFAULT)){
                  return true;
            }
	}
      return false;
}

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
  //the SnapshotName of illustris data is snapdir_XXX,which only has three letters between "snap" and "_"
  //the following lines add the eliminated "_" for illustris data
  if(subname[4]!='_'){
      subname.insert(4,"_");
  }
  stringstream formatter;
  formatter<<HBTConfig.SnapshotPath<<"/"<<SnapshotName<<"/"<<subname<<"."<<ifile<<".hdf5";
  filename=formatter.str();
}

void ApostleReader_t::GetIllustrisGroupName(int ifile, string &filename)
{
  string group_snap=SnapshotName;
  group_snap.replace(0,7,"groups"); //change snapdir_xxx to groups_xxx
  stringstream formatter_group;
  formatter_group<<HBTConfig.HaloPath<<"/"<<group_snap<<"/"<<group_snap<<"."<<ifile<<".hdf5";
  filename=formatter_group.str();
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
  is_illustris=IsIllustris(file);
  //PartType3 in illustris is TRACERS,which is not information of particles
  if(is_illustris){
      Header.npartTotal[3]=0;
  }
  H5Fclose(file);
}
void ApostleReader_t::GetParticleCountInFile(hid_t file, int np[])
{
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, np);
#ifdef DM_ONLY
  for(int i=0;i<TypeMax;i++)
	if(i!=TypeDM) np[i]=0;
#endif
  //PartType3 in illustris is TRACERS,which is not information of particles
  if(is_illustris){
      np[3]=0;
  }
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
	GetParticleCountInFile(file, np_this);
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
  GetParticleCountInFile(file, np_this.data());
  CompileOffsets(np_this, offset_this);

  HBTReal vunit=sqrt(Header.ScaleFactor);
  HBTReal boxsize=Header.BoxSize;
  for(int itype=0;itype<TypeMax;itype++)
  {
	int np=np_this[itype];
	if(np==0) continue;
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
	  if(H5Lexists(particle_data, "Velocities", H5P_DEFAULT))
	    ReadDataset(particle_data, "Velocities", H5T_HBTReal, v.data());
	  else
	    ReadDataset(particle_data, "Velocity", H5T_HBTReal, v.data());
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

#ifndef DM_ONLY
	//mass
	if(Header.mass[itype]==0)
	{
	  vector <HBTReal> m(np);
	  if(H5Lexists(particle_data, "Masses", H5P_DEFAULT))
	    ReadDataset(particle_data, "Masses", H5T_HBTReal, m.data());
	  else
	    ReadDataset(particle_data, "Mass", H5T_HBTReal, m.data());
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].Mass=m[i];
	}
	else
	{
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].Mass=Header.mass[itype];
	}

	//internal energy
#ifdef HAS_THERMAL_ENERGY
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
#endif

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
  GetParticleCountInFile(file, np_this.data());
  CompileOffsets(np_this, offset_this);

  for(int itype=0;itype<TypeMax;itype++)
  {
	int np=np_this[itype];
	if(np==0) continue;
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
		ParticlesThisType[i].HostId=(id[i]<0?NullGroupId:id[i]);//negative means particles inside virial radius but outside fof.
	}

	H5Gclose(particle_data);
  }

  H5Fclose(file);
}
void ApostleReader_t::ReadIllustrisGroup(int ifile, int itype,int &offset,ParticleHost_t *ParticlesInFile, bool FlagReadParticleId)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  vector <int> np_this(TypeMax);
  vector <HBTInt> offset_this(TypeMax);
  GetParticleCountInFile(file, np_this.data());
  CompileOffsets(np_this, offset_this);

  int np=np_this[itype];
  if(np!=0 ){
      auto ParticlesThisType=ParticlesInFile+offset;
      stringstream grpname;
      grpname<<"PartType"<<itype;
      if(H5Lexists(file, grpname.str().c_str(), H5P_DEFAULT)) {
            hid_t particle_data=H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);

            if(FlagReadParticleId)
            {//id
                  check_id_size(particle_data);
                  vector <HBTInt> id(np);
                  ReadDataset(particle_data, "ParticleIDs", H5T_HBTInt, id.data());
                  for(int i=0;i<np;i++)
                        ParticlesThisType[i].ParticleId=id[i];
            }

            H5Gclose(particle_data);
      }
  }
  offset=offset+np;//offset for this data type by file
  H5Fclose(file);
}

void ApostleReader_t::LoadSnapshot(int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology)
{
  SetSnapshot(snapshotId);
  ReadHeader(0, Header);
  Cosmology.Set(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);
  Cosmology.ParticleMass=Header.mass[TypeDM];
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
  typedef array <HBTInt, 6> HBTgroup;

  SetSnapshot(snapshotId);
  ReadHeader(0, Header);
  HBTInt NumberOfParticles=CompileFileOffsets(Header.NumberOfFiles);
  vector <ParticleHost_t> Particles(NumberOfParticles);
  bool FlagReadId=!HBTConfig.GroupLoadedIndex;

if(is_illustris){
//Offsets by type for all files
  vector <HBTInt> offsets(TypeMax);
  for (int itype=0;itype<TypeMax;itype++)
      {
            offsets[itype+1]=offsets[itype]+Header.npartTotal[itype];
      }

  cout<<"reading group files: ";
//read snapshot first by type ,then by file
#pragma omp parallel for num_threads(HBTConfig.MaxConcurrentIO)
  for (int itype=0;itype<TypeMax;itype++)
  {
      int off=0;
      for(int iFile=0; iFile<Header.NumberOfFiles; iFile++)
      {
            ReadIllustrisGroup(iFile,itype,off, Particles.data()+offsets[itype],1);
            if(itype==5){
            cout<<iFile<<" "<<flush;
            }
      }
  }
  cout<<endl;

//read head of group file
  string filename;
  GetIllustrisGroupName(0, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HBTInt NumGroups=0;
  HBTInt Numfiles=0;
  ReadAttribute(file, "Header", "Ngroups_Total", H5T_NATIVE_INT, &NumGroups);
  ReadAttribute(file, "Header", "NumFiles", H5T_NATIVE_INT, &Numfiles);
  H5Fclose(file);

//read group file
  vector <HBTgroup> group_offsets(NumGroups);
  HBTInt j=0;//the index for group_offset
  HBTInt j_original=0;
  for(HBTInt i=0;i<Numfiles;i++)
  {
	string filename;
	GetIllustrisGroupName(i, filename);
      hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      HBTInt NumGroups_thisfile=0;
      ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_NATIVE_INT, &NumGroups_thisfile);

	stringstream grpname;
	grpname<<"Offsets";
	hid_t group_data=H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);

      vector <HBTgroup> group_offsets_this(NumGroups_thisfile);
      ReadDataset(group_data, "Group_SnapByType",H5T_HBTInt, group_offsets_this.data());
      j_original=j;
      HBTInt k=0;
      while (j<j_original+NumGroups_thisfile)
      {
            group_offsets[j]=group_offsets_this[k];
            j=j+1;
            k=k+1;
      }

	H5Gclose(group_data);

      H5Fclose(file);
  }

  vector <HBTInt> group_offsets_alltype(NumGroups);
  HBTInt k=0;
  while (k<NumGroups)
      {
            group_offsets_alltype[k]=group_offsets[k][0]+group_offsets[k][1]+group_offsets[k][2]+group_offsets[k][3]+group_offsets[k][4]+group_offsets[k][5];
            k=k+1;
      }

  Halos.resize(NumGroups-1);
	for(HBTInt i=0;i<NumGroups-1;i++)
	{
	  HBTInt np=group_offsets_alltype[i+1]-group_offsets_alltype[i];
	  auto &p=Halos[i].Particles;
	  p.resize(np);
	  auto p_in=Particles.data();
	  //Offsets for single group by type
	  vector <HBTInt> single_group_offsets(TypeMax);
	  single_group_offsets[0]=group_offsets[i+1][0]-group_offsets[i][0];
	  for (int itype=1;itype<TypeMax;itype++)
        {
            single_group_offsets[itype]=single_group_offsets[itype-1]+group_offsets[i+1][itype]-group_offsets[i][itype];
        }
        auto p_in_o=p_in;
	  p_in=p_in+group_offsets[i][0];
	  int k=0;
	  for(HBTInt j=0;j<single_group_offsets[0];j++){
                  p[j]=p_in[k].ParticleId;
                  k=k+1;
            }
        for (int itype=1;itype<TypeMax;itype++)
        {
            p_in=p_in_o+offsets[itype]+group_offsets[i][itype];
            int k=0;
	      for(HBTInt j=single_group_offsets[itype-1];j<single_group_offsets[itype];j++){
                  p[j]=p_in[k].ParticleId;
                  k=k+1;
            }
	  }
	}

  cout<<Halos.size()<<" groups loaded";
  if(Halos.size()) cout<<" : "<<Halos[0].Particles.size();
  if(Halos.size()>1) cout<<","<<Halos[1].Particles.size()<<"...";
  cout<<endl;

  return NumberOfParticles;
}
else{
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

  MYSORT(Particles.begin(), Particles.end(), CompParticleHost);
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
}

bool IsApostleGroup(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 7)=="apostle";
}

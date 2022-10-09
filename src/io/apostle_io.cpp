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

#define PartTypeTracer 3

namespace Apostle
{
void ApostleSnap_t::SetSnapshot(int snapshotId)
{
  if(IsIllustrisSnap(HBTConfig.SnapshotFormat))
      SnapType=SnapshotType_t::illustris;
  if(IsApostleSnap(HBTConfig.SnapshotFormat))
      SnapType=SnapshotType_t::apostle;
  if(IsHelucidSnap(HBTConfig.SnapshotFormat))
      SnapType=SnapshotType_t::helucid;
  
//   SnapIsIllustris=IsIllustrisSnap(HBTConfig.SnapshotFormat);
  if(SnapType==SnapshotType_t::illustris||SnapType==SnapshotType_t::helucid)
    SnapDirBaseName="snapdir";
  else
    SnapDirBaseName="snapshot";
  
//   SnapshotId=snapshotId;
  if(HBTConfig.SnapshotNameList.empty())
  {
	stringstream formatter;
	formatter<<SnapDirBaseName<<"_"<<setw(3)<<setfill('0')<<snapshotId;
	SnapshotName=formatter.str();
  }
  else
	SnapshotName=HBTConfig.SnapshotNameList[snapshotId];
}

void ApostleHeader_t::GetFileName(int ifile, string &filename)
{
  string subname=SnapshotName;
  auto baselen=subname.find("_");
  subname.erase(4, baselen-4);//i.e., remove "shot" from "snapshot" or "snipshot"; or remove "dir" from "snapdir"
//   if(SnapType==SnapshotType_t::illustris||SnapType==SnapshotType_t::helucid)
//       subname.erase(4,3);//remove "dir" from "snapdir"
//   else
//       subname.erase(4, 4);//i.e., remove "shot" from "snapshot" or "snipshot"
  stringstream formatter;
  formatter<<HBTConfig.SnapshotPath<<"/"<<SnapshotName<<"/"<<subname<<"."<<ifile<<".hdf5";
  filename=formatter.str();
}

void ApostleHeader_t::GetParticleCountInFile(hid_t file, HBTInt np[])
{
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_HBTInt, np);
  //PartType3 in illustris is TRACERS,which is not information of particles
  if(SnapType==SnapshotType_t::illustris) np[PartTypeTracer]=0;
#ifdef DM_ONLY
  for(int i=0;i<TypeMax;i++)
	if(i!=TypeDM) np[i]=0;
#endif
}

void ApostleHeader_t::ReadFileHeader(int ifile)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  ReadAttribute(file, "Header", "NumFilesPerSnapshot", H5T_NATIVE_INT, &NumberOfFiles);
  ReadAttribute(file, "Header", "BoxSize", H5T_NATIVE_DOUBLE, &BoxSize);
  assert((HBTReal)BoxSize==HBTConfig.BoxSize);
  ReadAttribute(file, "Header", "Time", H5T_NATIVE_DOUBLE, &ScaleFactor);
  ReadAttribute(file, "Header", "Omega0", H5T_NATIVE_DOUBLE, &OmegaM0);
  ReadAttribute(file, "Header", "OmegaLambda", H5T_NATIVE_DOUBLE, &OmegaLambda0);
  ReadAttribute(file, "Header", "MassTable", H5T_NATIVE_DOUBLE, Mass);
//   cout<<Mass[0]<<","<<Mass[1]<<","<<Mass[2]<<endl;
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_HBTInt, NumPart.data());

  unsigned long long np[TypeMax], np_high[TypeMax];
  ReadAttribute(file, "Header", "NumPart_Total", H5T_NATIVE_ULLONG, np);
  ReadAttribute(file, "Header", "NumPart_Total_HighWord", H5T_NATIVE_ULLONG, np_high);
  for(int i=0;i<TypeMax;i++)
	NumPartTotal[i]=(np_high[i]<<32)|np[i];
  //Skip Tracers in Illustris. TODO: the gas particles in Illustris are moving cells, not completely Lagrangian. may need special treatment
  if(SnapType==SnapshotType_t::illustris) 
  {
    NumPart[PartTypeTracer]=0;
    NumPartTotal[PartTypeTracer]=0;
  }
  #ifdef DM_ONLY
  for(int i=0;i<TypeMax;i++)
	if(i!=TypeDM) NumPartTotal[i]=0;
  #endif
  H5Fclose(file);
  
  
  if(IsHelucidGroup(HBTConfig.GroupFileFormat))
  {
      int nbit=sizeof(HBTInt)*8;
      NullGroupId=((HBTInt)1)<<(nbit-2); //a large enough arbitrary number
      MinGroupId=1;
  }
  else//default to apostle; not used by illustris
  {
      NullGroupId=1<<30; //1073741824
      MinGroupId=0;
  }
  
}

/*
void ApostleHeader_t::LoadExtraHeaderParams()
{
    string filename=HBTConfig.ConfigFile;
    ifstream ifs;
    ifs.open(filename);
    if(!ifs.is_open())
    {
        cerr<<"Error opening parameter file "<<filename<<endl;
        exit(1);
    }
    string line;
    while(getline(ifs,line))
    {
        if(line.compare(0,13,"[ReaderExtra]")==0)//find the section
            break;
    }
    
    string name;
//     HBTInt NullGroupId, MinGroupId;
    #define ReadNameValue(x) ifs>>name>>x;assert(name==#x)
    ReadNameValue(NullGroupId);
    ReadNameValue(MinGroupId);
    #undef ReadNameValue
    ifs.close();
}
*/

void ApostleHeader_t::CompileTypeOffsets()
{
  NumPartAll=CompileOffsets(NumPartTotal.begin(), NumPartTotal.end(), TypeOffsetInMem.begin());
  
  NumPartEachFile.resize(NumberOfFiles);
  FileOffsetInType.resize(NumberOfFiles);
  
  TypeCounts_t offset_in_type={};
  
  for(int ifile=0;ifile<NumberOfFiles;ifile++)
  {
	auto &np_type=NumPartEachFile[ifile];
	string filename;
	GetFileName(ifile, filename);
	hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	GetParticleCountInFile(file, np_type.data());
	H5Fclose(file);
        
    FileOffsetInType[ifile]=offset_in_type;
    for(int itype=0;itype<TypeMax;itype++)
        offset_in_type[itype]+=np_type[itype];
  }
}

void ApostleHeader_t::Fill(int snapshotId)
{
    SetSnapshot(snapshotId);
    ReadFileHeader(0);
    CompileTypeOffsets();
}

void IllustrisGroupHeader_t::GetFileName(int ifile, string &filename)
{
  string snap_idname=SnapshotName.substr(SnapshotName.size()-3); //last 3 chars
  string group_snap="groups_"+snap_idname; 
  stringstream formatter, formatter2;
  formatter<<HBTConfig.HaloPath<<"/"<<group_snap<<"/fof_subhalo_tab_"<<snap_idname<<"."<<ifile<<".hdf5";
  formatter2<<HBTConfig.HaloPath<<"/"<<group_snap<<"/"<<group_snap<<"."<<ifile<<".hdf5"; //old illustris
  filename=formatter.str();
  if(!file_exist(filename.c_str()))
      filename=formatter2.str();
}

void IllustrisGroupHeader_t::ReadFileHeader(int ifile)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  ReadAttribute(file, "Header", "NumFiles", H5T_NATIVE_INT, &NumFiles);
  ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_HBTInt, &NumGroupsThisFile);
  ReadAttribute(file, "Header", "Ngroups_Total", H5T_HBTInt, &NumGroupsTotal);
  H5Fclose(file);
}

void IllustrisGroupHeader_t::CompileFileOffsets()
{
    GroupFileLen.resize(NumFiles);
    GroupFileOffset.resize(NumFiles);
    
    HBTInt ngroups, offset=0;
    for(int ifile=0;ifile<NumFiles;ifile++)
    {
        string filename;
        GetFileName(ifile, filename);
        hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_HBTInt, &ngroups);
        H5Fclose(file);
        
        GroupFileLen[ifile]=ngroups;
        GroupFileOffset[ifile]=offset;
        offset+=ngroups;
    }
    
    assert(offset==NumGroupsTotal);
}

void IllustrisGroupHeader_t::Fill(int snapshotId)
{
    SetSnapshot(snapshotId);
    ReadFileHeader(0);
    CompileFileOffsets();
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
void ApostleReader_t::ReadSnapshot(int ifile, Particle_t *Particles)
{
  string filename;
  SnapHeader.GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  auto &np_this=SnapHeader.NumPartEachFile[ifile];

  HBTReal vunit=sqrt(SnapHeader.ScaleFactor);
  HBTReal boxsize=SnapHeader.BoxSize;
  for(int itype=0;itype<TypeMax;itype++)
  {
	int np=np_this[itype];
	if(np==0) continue;
	auto ParticlesThisType=Particles+SnapHeader.TypeOffsetInMem[itype]+SnapHeader.FileOffsetInType[ifile][itype];//order type by type globally, and file by file in each type
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

#ifdef DM_ONLY
    //mass
    if(SnapHeader.Mass[TypeDM]==0)//mass not recorded in header
    {
        vector <HBTReal> m(np);
        if(H5Lexists(particle_data, "Masses", H5P_DEFAULT))
            ReadDataset(particle_data, "Masses", H5T_HBTReal, m.data());
        else
            ReadDataset(particle_data, "Mass", H5T_HBTReal, m.data());
        assert(m.front()==m.back()); //must be uniform mass
#pragma omp single nowait
        SnapHeader.Mass[TypeDM]=m.front();
        cout<<"WARNING: running in DM_ONLY mode, but header does not contain a fixed DM particle mass.\n Setting to first mass = "<<m.front()<<" from mass record"<<endl;
    }
#else
	//mass
	if(SnapHeader.Mass[itype]==0)
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
		ParticlesThisType[i].Mass=SnapHeader.Mass[itype];
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

void ApostleReader_t::ReadGroupId(int ifile, ParticleHost_t *Particles, bool FlagReadParticleId)
{//read apostle group HostIds
  string filename;
  SnapHeader.GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  auto &np_this=SnapHeader.NumPartEachFile[ifile];

  for(int itype=0;itype<TypeMax;itype++)
  {
	int np=np_this[itype];
	if(np==0) continue;
	auto ParticlesThisType=Particles+SnapHeader.TypeOffsetInMem[itype]+SnapHeader.FileOffsetInType[ifile][itype];//order type by type globally, and file by file in each type
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
      if(H5Lexists(particle_data, "GroupNumber", H5P_DEFAULT))
          ReadDataset(particle_data, "GroupNumber", H5T_HBTInt, id.data());
      else
          ReadDataset(particle_data, "HaloID", H5T_HBTInt, id.data());
	  for(int i=0;i<np;i++)
      {
          HBTInt id_fixed=id[i]-SnapHeader.MinGroupId;//shift the HaloId to start from 0 if not.
          ParticlesThisType[i].HostId=(id_fixed<0?SnapHeader.NullGroupId:id_fixed);//negative means particles inside virial radius but outside fof.
      }
	}

	H5Gclose(particle_data);
  }

  H5Fclose(file);
}

void ApostleReader_t::ReadParticleIDs(int ifile, HBTInt *ParticleIDs)
{//Read ID only
  string filename;
  SnapHeader.GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  auto &np_this=SnapHeader.NumPartEachFile[ifile];

  for(int itype=0;itype<TypeMax;itype++)
  {
	int np=np_this[itype];
	if(np==0) continue;
	auto ParticlesThisType=ParticleIDs+SnapHeader.TypeOffsetInMem[itype]+SnapHeader.FileOffsetInType[ifile][itype];//order type by type globally, and file by file in each type
	stringstream grpname;
	grpname<<"PartType"<<itype;
	if(!H5Lexists(file, grpname.str().c_str(), H5P_DEFAULT)) continue;
	hid_t particle_data=H5Gopen2(file, grpname.str().c_str(), H5P_DEFAULT);
// 	if(particle_data<0) continue; //skip non-existing type

	check_id_size(particle_data);
    ReadDataset(particle_data, "ParticleIDs", H5T_HBTInt, ParticlesThisType);

	H5Gclose(particle_data);
  }

  H5Fclose(file);
}

void ApostleReader_t::LoadSnapshot(int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology)
{
  SnapHeader.Fill(snapshotId);
  Cosmology.Set(SnapHeader.ScaleFactor, SnapHeader.OmegaM0, SnapHeader.OmegaLambda0);
  
  Particles.resize(SnapHeader.NumPartAll);
#pragma omp parallel for num_threads(HBTConfig.MaxConcurrentIO)
  for(int iFile=0; iFile<SnapHeader.NumberOfFiles; iFile++)
  {
	ReadSnapshot(iFile, Particles.data());
	cout<<iFile<<" "<<flush;
  }
  cout<<endl;
  cout<<" ( "<<SnapHeader.NumberOfFiles<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
  
  Cosmology.ParticleMass=SnapHeader.Mass[TypeDM];
//   cout<<" Particle[0]: x="<<Particles[0].ComovingPosition<<", v="<<Particles[0].PhysicalVelocity<<", m="<<Particles[0].Mass<<endl;
//   cout<<" Particle[2]: x="<<Particles[2].ComovingPosition<<", v="<<Particles[2].PhysicalVelocity<<", m="<<Particles[2].Mass<<endl;
//   cout<<" Particle[end]: x="<<Particles.back().ComovingPosition<<", v="<<Particles.back().PhysicalVelocity<<", m="<<Particles.back().Mass<<endl;
}

inline bool CompParticleHost(const ParticleHost_t &a, const ParticleHost_t &b)
{
  return a.HostId<b.HostId;
}

HBTInt ApostleReader_t::LoadApostleGroups(int snapshotId, vector< Halo_t >& Halos)
{
  SnapHeader.Fill(snapshotId);
  vector <ParticleHost_t> Particles(SnapHeader.NumPartAll);
  bool FlagReadId=!HBTConfig.GroupLoadedIndex;

  cout<<"reading group files: ";
  #pragma omp parallel for num_threads(HBTConfig.MaxConcurrentIO)
  for(int iFile=0; iFile<SnapHeader.NumberOfFiles; iFile++)
  {
	ReadGroupId(iFile, Particles.data(), FlagReadId);
	cout<<iFile<<" "<<flush;
  }
  cout<<endl;

  if(!FlagReadId)
  {//fill with index
	#pragma omp parallel for
	for(HBTInt i=0;i<SnapHeader.NumPartAll;i++)
	  Particles[i].ParticleId=i;
  }

  MYSORT(Particles.begin(), Particles.end(), CompParticleHost);
  assert(Particles.back().HostId==SnapHeader.NullGroupId);//max haloid==NullGroupId
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

HBTInt ApostleReader_t::LoadIllustrisGroups(int snapshotId, vector< Halo_t >& Halos)
{
//   typedef array <HBTInt, 6> HBTgroup;

  SnapHeader.Fill(snapshotId);
  GroupHeader.Fill(snapshotId);
  Halos.resize(GroupHeader.NumGroupsTotal);
  
  vector <TypeCounts_t> GroupLenType, GroupOffsetType;
  GroupLenType.resize(GroupHeader.NumGroupsTotal);
  GroupOffsetType.resize(GroupHeader.NumGroupsTotal);
  
  #pragma omp parallel for num_threads(HBTConfig.MaxConcurrentIO)
  for(int ifile=0;ifile<GroupHeader.NumFiles;ifile++)
  {
    if(GroupHeader.GroupFileLen[ifile]==0) continue; //skip empty files
    string filename;
    GroupHeader.GetFileName(ifile, filename);
    hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ReadDataset(file, "Group/GroupLenType", H5T_HBTInt, GroupLenType.data()+GroupHeader.GroupFileOffset[ifile]);
    H5Fclose(file);
   }
   
   
  #pragma omp parallel for
  for(int igrp=0;igrp<GroupHeader.NumGroupsTotal;igrp++)
  {
    GroupLenType[igrp][PartTypeTracer]=0;//skip tracers
    int np=0;
    for(int itype=0;itype<TypeMax;itype++)
    {
        #ifdef DM_ONLY
        if(itype!=TypeDM) GroupLenType[igrp][itype]=0; //clear other types
        #endif
        np+=GroupLenType[igrp][itype];
    }
    Halos[igrp].Particles.resize(np);
  }
  
  //TODO: test if offset file exist, if yes, read it to fill the offset
  
  TypeCounts_t type_offset=SnapHeader.TypeOffsetInMem;
  #pragma omp parallel for
  for(int itype=0;itype<TypeMax;itype++)
  {
    for(int igrp=0;igrp<GroupHeader.NumGroupsTotal;igrp++)
    {
        GroupOffsetType[igrp][itype]=type_offset[itype];
        type_offset[itype]+=GroupLenType[igrp][itype];
    }
  }
  
  #pragma omp parallel for
  for(int igrp=0;igrp<GroupHeader.NumGroupsTotal;igrp++)
  {
    int offset_local=0;
    for(int itype=0;itype<TypeMax;itype++)
    {
        int np_thistype=GroupLenType[igrp][itype];
        HBTInt pid=GroupOffsetType[igrp][itype];
        auto particles=Halos[igrp].Particles.data()+offset_local;
        offset_local+=np_thistype;
        for(int i=0;i<np_thistype;i++)
            particles[i]=pid++;
    }
  }
    
  
  bool FlagReadId=!HBTConfig.GroupLoadedIndex;
  if(FlagReadId)
  {
    vector <HBTInt> ParticleIDs(SnapHeader.NumPartAll);
    #pragma omp parallel for num_threads(HBTConfig.MaxConcurrentIO)
    for(int iFile=0; iFile<SnapHeader.NumberOfFiles; iFile++)
        ReadParticleIDs(iFile, ParticleIDs.data());
        
    #pragma omp parallel for 
    for(int igrp=0; igrp<Halos.size();igrp++)
    {
        auto &particles=Halos[igrp].Particles;
        for(int i=0;i<particles.size();i++)
            particles[i]=ParticleIDs[particles[i]]; //convert index to id
    }
  }

  cout<<Halos.size()<<" groups loaded";
  if(Halos.size()) cout<<", np= "<<Halos[0].Particles.size();
  if(Halos.size()>1) cout<<","<<Halos[1].Particles.size()<<"...";
  cout<<endl;

  return Halos.size();
}


bool IsHelucidGroup(const string& GroupFileFormat)
{//helucid is a minor variant of apostle
    return GroupFileFormat.substr(0,15)=="apostle_helucid";
}
bool IsHelucidSnap(const string &SnapshotFormat)
{
    return SnapshotFormat.substr(0,15)=="apostle_helucid";
}
bool IsApostleGroup(const string &GroupFileFormat)
{//applies to both apostle and helucid; 
    return GroupFileFormat.substr(0, 7)=="apostle";
}
bool IsIllustrisGroup(const string &GroupFileFormat)
{
    return GroupFileFormat.substr(0,9)=="illustris";
}
bool IsApostleSnap(const string &SnapshotFormat)
{//applies to both apostle and helucid
    return SnapshotFormat.substr(0,7)=="apostle";
}
bool IsIllustrisSnap(const string &SnapshotFormat)
{
    return SnapshotFormat.substr(0,9)=="illustris";
}
}

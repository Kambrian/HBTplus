using namespace std;
#include <iostream>
#include <numeric>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <list>
#include <stdexcept>

#include "../snapshot.h"
#include "../mymath.h"
#include "../hdf_wrapper.h"
#include "gadget4_io.h"

namespace Gadget4Reader
{
void create_Gadget4Header_MPI_type(MPI_Datatype& dtype)
{
  /*to create the struct data type for communication*/
  Gadget4Header_t p;
  #define NumAttr 13
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint   offsets[NumAttr], origin,extent;

  MPI_Get_address(&p,&origin);
  MPI_Get_address((&p)+1,&extent);//to get the extent of s
  extent-=origin;

  int i=0;
  #define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
  RegisterAttr(NumberOfFiles, MPI_INT, 1)
  RegisterAttr(BoxSize, MPI_DOUBLE, 1)
  RegisterAttr(ScaleFactor, MPI_DOUBLE, 1)
  RegisterAttr(BoxSize, MPI_DOUBLE, 1)
  RegisterAttr(OmegaM0, MPI_DOUBLE, 1)
  RegisterAttr(OmegaLambda0, MPI_DOUBLE, 1)
  RegisterAttr(mass, MPI_DOUBLE, TypeMax)
  RegisterAttr(npart[0], MPI_INT, TypeMax)
  RegisterAttr(npartTotal[0], MPI_HBT_INT, TypeMax)
  #undef RegisterAttr
  assert(i<=NumAttr);

  MPI_Type_create_struct(i,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef NumAttr
}

void Gadget4Reader_t::SetSnapshot(int snapshotId)
{
  if(HBTConfig.SnapshotNameList.empty())
  {
	stringstream formatter;
	formatter<<"snapdir_"<<setw(3)<<setfill('0')<<snapshotId;
	SnapshotName=formatter.str();
  }
  else
	SnapshotName=HBTConfig.SnapshotNameList[snapshotId];
}

void Gadget4Reader_t::GetFileName(int ifile, string &filename)
{
  string snap_idname=SnapshotName.substr(SnapshotName.size()-3); //last 3 chars
  stringstream formatter;
  formatter<<HBTConfig.SnapshotPath<<"/"<<SnapshotName<<"/"<<HBTConfig.SnapshotFileBase<<"_"<<snap_idname<<"."<<ifile<<".hdf5";
  filename=formatter.str();
}

void Gadget4Reader_t::GetGroupFileName(int ifile, string &filename)
{
  string snap_idname=SnapshotName.substr(SnapshotName.size()-3); //last 3 chars
  stringstream formatter;
  formatter<<HBTConfig.HaloPath<<"/groups_"<<snap_idname<<"/fof_subhalo_tab_"<<snap_idname<<"."<<ifile<<".hdf5";
  filename=formatter.str();
}

void Gadget4Reader_t::ReadHeader(int ifile, Gadget4Header_t &header)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  ReadAttribute(file, "Header", "NumFilesPerSnapshot", H5T_NATIVE_INT, &header.NumberOfFiles);
  ReadAttribute(file, "Header", "BoxSize", H5T_NATIVE_DOUBLE, &header.BoxSize);
  assert((HBTReal)header.BoxSize==HBTConfig.BoxSize);
  ReadAttribute(file, "Header", "Time", H5T_NATIVE_DOUBLE, &header.ScaleFactor);
  ReadAttribute(file, "Header", "Omega0", H5T_NATIVE_DOUBLE, &header.OmegaM0);
  ReadAttribute(file, "Header", "OmegaLambda", H5T_NATIVE_DOUBLE, &header.OmegaLambda0);
  ReadAttribute(file, "Header", "MassTable", H5T_NATIVE_DOUBLE, header.mass);
  cout<<"mass table: "<<header.mass[0]<<","<<header.mass[1]<<","<<header.mass[2]<<endl;
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, header.npart);

//   unsigned np[TypeMax], np_high[TypeMax];
  ReadAttribute(file, "Header", "NumPart_Total", H5T_HBTInt, header.npartTotal);
//   ReadAttribute(file, "Header", "NumPart_Total_HighWord", H5T_NATIVE_UINT, np_high);
//   for(int i=0;i<TypeMax;i++)
// 	header.npartTotal[i]=(((unsigned long)np_high[i])<<32)|np[i];
  H5Fclose(file);
}

int Gadget4Reader_t::ReadGroupFileCounts(int ifile)
{
  int nfiles;
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  ReadAttribute(file, "Header", "NumFiles", H5T_NATIVE_INT, &nfiles);
  H5Fclose(file);

  return nfiles;
}

void Gadget4Reader_t::GetParticleCountInFile(hid_t file, int np[])
{
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, np);
#ifdef DM_ONLY
  for(int i=0;i<TypeMax;i++)
	if(i!=TypeDM) np[i]=0;
#endif
}
HBTInt Gadget4Reader_t::CompileFileOffsets(int nfiles)
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

HBTInt Gadget4Reader_t::CompileGroupFileOffsets(int nfiles)
{
  HBTInt offset=0, nhalo_total;
  nhalo_per_groupfile.reserve(nfiles);
  offsethalo_per_groupfile.reserve(nfiles);
  for(int ifile=0;ifile<nfiles;ifile++)
  {
	offsethalo_per_groupfile.push_back(offset);

	HBTInt nhalo_this;
	string filename;
	GetGroupFileName(ifile, filename);
	hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_HBTInt, &nhalo_this);
    if(ifile==0)
      ReadAttribute(file, "Header", "Ngroups_Total", H5T_HBTInt, &nhalo_total);
	H5Fclose(file);

	nhalo_per_groupfile.push_back(nhalo_this);
	offset+=nhalo_this;
  }
  assert(nhalo_total==offset);
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
void Gadget4Reader_t::ReadSnapshot(int ifile, Particle_t *ParticlesInFile)
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

#ifndef DM_ONLY
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

void Gadget4Reader_t::ReadGroupLen(int ifile, HBTInt *buf)
{
  string filename;
  GetFileName(ifile, filename);
  hid_t file=H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HBTInt ngroups;
  ReadAttribute(file, "Header", "Ngroups_ThisFile", H5T_HBTInt, &ngroups);
  if(ngroups>0)
    ReadDataset(file, "Group/GroupLen", H5T_HBTInt, buf);
  H5Fclose(file);
}

void Gadget4Reader_t::ReadGroupParticles(int ifile, ParticleHost_t *ParticlesInFile, bool FlagReadParticleId)
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

	if(FlagReadParticleId)
	{
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

#ifndef DM_ONLY
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
	}

	{//Hostid
	  vector <HBTInt> id(np);
	  ReadDataset(particle_data, "GroupNumber", H5T_HBTInt, id.data());
	  for(int i=0;i<np;i++)
		ParticlesThisType[i].HostId=(id[i]<0?NullGroupId:id[i]);//negative means outside fof but within Rv
	}

	H5Gclose(particle_data);
  }

  H5Fclose(file);
}

void Gadget4Reader_t::LoadSnapshot(MpiWorker_t &world, int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology)
{
  SetSnapshot(snapshotId);

  const int root=0;
  if(world.rank()==root)
  {
    ReadHeader(0, Header);
    CompileFileOffsets(Header.NumberOfFiles);
  }
  MPI_Bcast(&Header, 1, MPI_Gadget4Header_t, root, world.Communicator);
  world.SyncContainer(np_file, MPI_HBT_INT, root);
  world.SyncContainer(offset_file, MPI_HBT_INT, root);

  Cosmology.Set(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);
#ifdef DM_ONLY
//   Cosmology.ParticleMass=Header.mass[TypeDM];
#endif

  HBTInt nfiles_skip, nfiles_end;
  AssignTasks(world.rank(), world.size(), Header.NumberOfFiles, nfiles_skip, nfiles_end);
  {
    HBTInt np=0;
    np=accumulate(np_file.begin()+nfiles_skip, np_file.begin()+nfiles_end, np);
    Particles.resize(np);
  }

  for(int i=0, ireader=0;i<world.size();i++, ireader++)
  {
	if(ireader==HBTConfig.MaxConcurrentIO)
	{
	  ireader=0;//reset reader count
	  MPI_Barrier(world.Communicator);//wait for every thread to arrive.
	}
	if(i==world.rank())//read
	{
	  for(int iFile=nfiles_skip; iFile<nfiles_end; iFile++)
	  {
// 	        cout<<"("<<world.rank()<<","<<iFile<<")"<<flush;
		ReadSnapshot(iFile, Particles.data()+offset_file[iFile]-offset_file[nfiles_skip]);
	  }
	}
  }

//   cout<<endl;
//   cout<<" ( "<<Header.NumberOfFiles<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
//   cout<<" Particle[0]: x="<<Particles[0].ComovingPosition<<", v="<<Particles[0].PhysicalVelocity<<", m="<<Particles[0].Mass<<endl;
//   cout<<" Particle[2]: x="<<Particles[2].ComovingPosition<<", v="<<Particles[2].PhysicalVelocity<<", m="<<Particles[2].Mass<<endl;
//   cout<<" Particle[end]: x="<<Particles.back().ComovingPosition<<", v="<<Particles.back().PhysicalVelocity<<", m="<<Particles.back().Mass<<endl;
}

typedef vector <HBTInt> CountBuffer_t;

struct FileAssignment_t
{
  int ifile_begin, ifile_end;
  HBTInt firstfile_begin_part, lastfile_end_part, npart;
  HBTInt haloid_begin, nhalo;
  FileAssignment_t()
  {
	ifile_begin=0;
	ifile_end=0;
	firstfile_begin_part=0;
	lastfile_end_part=-1;
	npart=0;
	haloid_begin=0;
	nhalo=0;
  }
  void create_MPI_type(MPI_Datatype &dtype)
  /*to create the struct data type for communication*/
  {
	FileAssignment_t &p=*this;
	#define NumAttr 7
	MPI_Datatype oldtypes[NumAttr];
	int blockcounts[NumAttr];
	MPI_Aint   offsets[NumAttr], origin,extent;

	MPI_Get_address(&p,&origin);
	MPI_Get_address((&p)+1,&extent);//to get the extent of s
	extent-=origin;

	int i=0;
	#define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
	RegisterAttr(ifile_begin, MPI_INT, 1)
	RegisterAttr(ifile_end, MPI_INT, 1)
	RegisterAttr(firstfile_begin_part, MPI_HBT_INT, 1)
	RegisterAttr(lastfile_end_part, MPI_HBT_INT, 1)
	RegisterAttr(npart, MPI_HBT_INT, 1)
	RegisterAttr(haloid_begin, MPI_HBT_INT, 1)
	RegisterAttr(nhalo, MPI_HBT_INT, 1)
	#undef RegisterAttr
	assert(i==NumAttr);

	MPI_Type_create_struct(NumAttr,blockcounts,offsets,oldtypes, &dtype);
	MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
	MPI_Type_commit(&dtype);
	#undef NumAttr
  }
};
typedef vector <HBTInt> ParticleIdBuffer_t;
typedef vector <HBTInt> CountBuffer_t;
void AssignHaloTasks(int nworkers, HBTInt npart_tot, const CountBuffer_t & HaloOffsets, const CountBuffer_t &FileOffsets, HBTInt & npart_begin, FileAssignment_t &task)
{
  if(0==npart_tot) return;

  HBTInt npart_this=(npart_tot-npart_begin)/nworkers;
  HBTInt npart_end=npart_begin+npart_this;

  task.haloid_begin=lower_bound(HaloOffsets.begin(), HaloOffsets.end()-1, npart_begin)-HaloOffsets.begin();
  HBTInt endhalo=lower_bound(HaloOffsets.begin(), HaloOffsets.end()-1, npart_end)-HaloOffsets.begin();
  task.nhalo=endhalo-task.haloid_begin;

  npart_end=HaloOffsets.at(endhalo);//need to append np to HaloOffsets to avoid overflow.
  task.npart=npart_end-npart_begin;

  task.ifile_begin=upper_bound(FileOffsets.begin(), FileOffsets.end(), npart_begin)-1-FileOffsets.begin();
  task.ifile_end=lower_bound(FileOffsets.begin(), FileOffsets.end(), npart_end)-FileOffsets.begin();

  task.firstfile_begin_part=npart_begin-FileOffsets[task.ifile_begin];
  task.lastfile_end_part=npart_end-FileOffsets[task.ifile_end-1];

  npart_begin=npart_end;
}

struct HaloProcMapper_t
{//metadata to facilitate the reading of local halo particles from the local snapshot on each process
    HBTInt Nhalos=0;
    HBTInt UpcomingHaloID=0; //the first halo located to the right of the current proc
    HBTInt UpcomingHaloPartOffset=0;//local offset of upcoming halo
    HBTInt PrevHaloLocalLen=0;//number of particles to send to previous halo
    int PrevHaloProcId=0; //processor on which the previous halo first appeared, so that we know where to send data
    vector<HBTInt> LastHaloProcLen;//segment sizes of the final halo on different procs
    vector<HBTInt> HaloLen;

    void FillUpcomingHalo(HBTInt haloid, HBTInt offset)
    {
      UpcomingHaloID=haloid;
      UpcomingHaloPartOffset=offset;
    }
    void FillRemaining(int NextProcId, HBTInt NextUpcomingHaloID, const vector<HBTInt> &HaloLenAll, const vector<HBTInt> &HaloOffsetAll, const vector <HBTInt> &ProcOffset)
    {
      Nhalos=NextUpcomingHaloID-UpcomingHaloID;//will be 0 if no halo actually starts on this proc
      HaloLen.assign(HaloLenAll.begin()+UpcomingHaloID, HaloLenAll.begin()+NextUpcomingHaloID);
      assert(LastHaloProcLen.size()==0); //make sure it is default initialized to empty
      if(Nhalos)
      {//fill segment sizes of LastHalo
        HBTInt offset_origin=HaloOffsetAll[NextUpcomingHaloID-1];
        for(int iproc_fast=NextProcId; iproc_fast<ProcOffset.size(); iproc_fast++)//find all upcoming segments
        {
          if(ProcOffset[iproc_fast]>=HaloOffsetAll[NextUpcomingHaloID])//past next halo
          {
            LastHaloProcLen.push_back(HaloOffsetAll[NextUpcomingHaloID]-offset_origin);//last proc
            break;
          }
          else
            LastHaloProcLen.push_back(ProcOffset[iproc_fast]-offset_origin);
          offset_origin=ProcOffset[iproc_fast];
        }
      }
    }
    void create_MPI_Shell_Type(MPI_Datatype &dtype)
  /*to create the struct data type for communication; vector members are not included.*/
  {
	FileAssignment_t &p=*this;
	#define NumAttr 5
	MPI_Datatype oldtypes[NumAttr];
	int blockcounts[NumAttr];
	MPI_Aint   offsets[NumAttr], origin,extent;

	MPI_Get_address(&p,&origin);
	MPI_Get_address((&p)+1,&extent);//to get the extent of s
	extent-=origin;

	int i=0;
	#define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
	RegisterAttr(Nhalos, MPI_HBT_INT, 1)
	RegisterAttr(UpcomingHaloID, MPI_HBT_INT, 1)
	RegisterAttr(UpcomingHaloPartOffset, MPI_HBT_INT, 1)
    RegisterAttr(PrevHaloLocalLen, MPI_HBT_INT, 1)
    RegisterAttr(PrevHaloProcId, MPI_INT, 1)
	#undef RegisterAttr
	assert(i==NumAttr);

	MPI_Type_create_struct(NumAttr,blockcounts,offsets,oldtypes, &dtype);
	MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
	MPI_Type_commit(&dtype);
	#undef NumAttr
  }
};

void Gadget4Reader_t::LoadGroups(MpiWorker_t &world, const ParticleSnapshot_t &partsnap, vector <Halo_t> &Halos)
{
  SetSnapshot(partsnap.GetSnapshotId());

  const auto &Particles=partsnap.Particles;

  int root=0;

  /* read file metadata */
  int FileCounts;
  HBTInt NhaloTotal;
  if(world.rank()==root)
  {
    FileCounts=ReadGroupFileCounts(0);
    NhaloTotal=CompileGroupFileOffsets(FileCounts);
  }
  world.SyncAtom(&FileCounts, MPI_INT, root);
  world.SyncAtom(&NhaloTotal, MPI_HBT_INT, root);
  world.SyncContainer(nhalo_per_groupfile, MPI_HBT_INT, root);
  world.SyncContainer(offsethalo_per_groupfile, MPI_HBT_INT, root);

  /* read halolen in parallel*/
  CountBuffer_t HaloLenLocal, HaloLenAll;
  HBTInt nfiles_skip, nfiles_end;
  AssignTasks(world.rank(), world.size(), FileCounts, nfiles_skip, nfiles_end);
  {
    int nhalo_local=0;
    nhalo_local=accumulate(nhalo_per_groupfile.begin()+nfiles_skip, nhalo_per_groupfile.begin()+nfiles_end, nhalo_local);
    HaloLenLocal.resize(nhalo_local);

    for(int i=0, ireader=0;i<world.size();i++, ireader++)
    {
      if(ireader==HBTConfig.MaxConcurrentIO)
      {
        ireader=0;//reset reader count
        MPI_Barrier(world.Communicator);//wait for every thread to arrive.
      }
      if(i==world.rank())//read
      {
        for(int iFile=nfiles_skip; iFile<nfiles_end; iFile++)
        {
          if(nhalo_per_groupfile[iFile])//some files do not have groups
            ReadGroupLen(iFile, HaloLenLocal.data()+offsethalo_per_groupfile[iFile]-offsethalo_per_groupfile[nfiles_skip]);
        }
      }
    }

    vector<int> buffersizes;
    if(world.rank()==root)
    {
      HaloLenAll.resize(NhaloTotal);
      buffersizes.resize(world.size());
    }
    MPI_Gather(&nhalo_local, 1, MPI_INT, buffersizes.data(), 1, MPI_INT, root, world.Communicator);
    MPI_Gatherv(HaloLenLocal.data(), HaloLenLocal.size(), MPI_HBT_INT, HaloLenAll.data(), buffersizes.data(), MPI_HBT_INT, root, world.Communicator);
  }

  //distribute read tasks
  CountBuffer_t ProcLen;
  {
    if(world.rank()==root) ProcLen.resize(world.size());
    HBTInt NumPartThisProc=Particles.size();
    MPI_Gather(&NumPartThisProc, 1, MPI_HBT_INT, ProcLen.data(), 1, MPI_HBT_INT, root, world.Communicator);
  }

  //build HaloMapper, for distributed reading of halo particles
  vector<HaloProcMapper_t> HaloMapperArr;
  if(world.rank()==root)
  {
    HaloMapperArr.resize(world.size());
    CountBuffer_t ProcOffset, HaloOffsetAll;
    HBTInt NumPartAllProc=CompileOffsets(ProcLen, ProcOffset);
    ProcOffset.push_back(NumPartAllProc);
    HBTInt NumPartAllHalos=CompileOffsets(HaloLenAll, HaloOffsetAll);
    HaloOffsetAll.push_back(NumPartAllHalos);
    cout<<"Fraction of Particles in halos: "<<NumPartAllHalos*1./NumPartAllProc;
    int iproc=0;
    for(HBTInt ihalo=0;ihalo<HaloOffsetAll.size();ihalo++)
    {
      while(HaloOffsetAll[ihalo]>=ProcOffset[iproc])//iterate over procs facing this upcoming halo
      {
        if(ProcLen.size()==iproc)
        {
          throw runtime_error("Offset of halo "+to_string(ihalo)+" is above snapshotsize "+to_string(ProcOffset[iproc])+endl)
          exit(1);
        };

        HaloMapperArr[iproc].FillUpcomingHalo(ihalo, HaloOffsetAll[ihalo]-ProcOffset[iproc]);
        //offset will exceed ProcLen if not a local halo
        if(iproc>0) //complete previous mapper
        {
          HaloMapperArr[iproc-1].FillRemaining(iproc, ihalo, HaloLenAll, HaloOffsetAll, ProcOffset);
          auto mapper=HaloMapperArr.begin()+iproc-1;//inform connected mappers
          for(int i=1;i<mapper->LastHaloProcLen.size();i++)
          {
            (mapper+i)->PrevHaloProcId=iproc-1;
            (mapper+i)->PrevHaloLocalLen=mapper->LastHaloProcLen[i];
          }
        }
        iproc++;
      }
    }
    if(iproc>0)
      HaloMapperArr[iproc-1].FillRemaining(iproc, HaloLenAll.size(), HaloLenAll, HaloOffsetAll, ProcOffset);
  }

  //distribute halo-mapper
  HaloProcMapper_t HaloMapperLocal;
  {
    MPI_Datatype MPI_HaloMapper_Shell_t;
    HaloProcMapper_t().create_MPI_Shell_Type(MPI_HaloMapper_Shell_t);
    MPI_Scatter(HaloMapperArr.data(), 1, MPI_HaloMapper_Shell_t, &HaloMapperLocal, 1, MPI_HaloMapper_Shell_t, root, world.Communicator);
    MPI_Type_free(&MPI_HaloMapper_Shell_t);
    int LastHaloSpan=0;
    vector <int> LastHaloSpanArr(world.size());
    for(int rank=0; rank<world.size(); rank++)
      LastHaloSpanArr[rank]=HaloMapperArr[rank].LastHaloProcLen.size();
    MPI_Scatter(LastHaloSpanArr.data(), 1, MPI_INT, &LastHaloSpan, 1, MPI_INT, root, world.Communicator);
    HaloMapperLocal.LastHaloProcLen.resize(LastHaloSpan);
    HaloMapperLocal.HaloLen.resize(HaloMapperLocal.Nhalos);
    if(world.rank()==root)
    {
      for(int rank=0;rank<world.size();rank++)
      {
        MPI_Send(HaloMapperArr[rank].LastHaloProcLen.data(), LastHaloSpanArr[rank], MPI_HBT_INT, rank, 0, world.Communicator);
        MPI_Send(HaloMapperArr[rank].HaloLen.data(), HaloMapperArr[rank].Nhalos, MPI_HBT_INT, rank, 1, world.Communicator);
      }
    }
    MPI_Recv(HaloMapperLocal.LastHaloProcLen.data(), LastHaloSpan, MPI_INT, root, 0, world.Communicator, MPI_STATUS_IGNORE);
    MPI_Recv(HaloMapperLocal.HaloLen.data(), HaloMapperLocal.Nhalos, MPI_INT, root, 1, world.Communicator, MPI_STATUS_IGNORE);
  }

  //read halo from processors
  Halos.resize(HaloMapperLocal.Nhalos);
  for(HBTInt i=0;i<Halos.size();i++)
  {
    Halos[i].HaloId=i+HaloMapperLocal.UpcomingHaloID;
    Halos[i].Particles.resize(HaloMapperLocal.HaloLen[i]);
  }
  //copy local particles to halo
  auto it_part=Particles.begin()+HaloMapperLocal.UpcomingHaloPartOffset;
  vector <HBTInt> HaloOffsetsLocal;
  CompileOffsets(HaloMapperLocal.HaloLen, HaloOffsetsLocal);
  for(HBTInt i=0;i<Halos.size()-1;i++)
    Halos[i].Particles.assign(it_part+HaloOffsetsLocal[i], it_part+HaloOffsetsLocal[i+1]);

//send out remote particles
  MPI_Request req_to_prev;
  if(HaloMapperLocal.PrevHaloLocalLen)
    MPI_Isend(Particles.data(), HaloMapperLocal.PrevHaloLocalLen, MPI_HBT_Particle, HaloMapperLocal.PrevHaloProcId, 2, world.Communicator, &req_to_prev);
//populate last halo
  auto &h=Halos.back();
  if(Halos.size()>0)
  {
    //populate local particles for last halo
    it_part+=HaloOffsetsLocal.back();
    assert(HaloMapperLocal.LastHaloProcLen.size()>0);
    for(HBTInt i=0;i<HaloMapperLocal.LastHaloProcLen[0];i++)
      h.Particles[i]=*(it_part+i);
    //fetch remote particles for last halo
    auto it_halopart=h.Particles.begin()+HaloMapperLocal.LastHaloProcLen[0];
    vector<HBTInt> LastHaloProcOffset;
    CompileOffsets(HaloMapperLocal.LastHaloProcLen, LastHaloProcOffset);
    for(int iproc=1;iproc<HaloMapperLocal.LastHaloProcLen.size();iproc++)
    {
      MPI_Recv(h.Particles.data()+LastHaloProcOffset[iproc], HaloMapperLocal.LastHaloProcLen[iproc], MPI_HBT_Particle, world.rank()+i, 2, world.Communicator, MPI_STATUS_IGNORE);
    }
  }
  MPI_Wait(&req_to_prev, MPI_STATUS_IGNORE);

  {//exchange halos to place them into subboxes
  MPI_Datatype MPI_HaloShell_t;
  create_MPI_Halo_Id_type(MPI_HaloShell_t);
  vector <Halo_t> LocalHalos;
  partsnap.ExchangeHalos(world, Halos, LocalHalos, MPI_HaloShell_t, false);
  Halos.swap(LocalHalos);
  My_Type_free(MPI_HaloShell_t);
  }
//   cout<<endl;
//   cout<<" ( "<<Header.NumberOfFiles<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
//   cout<<" Particle[0]: x="<<Particles[0].ComovingPosition<<", v="<<Particles[0].PhysicalVelocity<<", m="<<Particles[0].Mass<<endl;
//   cout<<" Particle[2]: x="<<Particles[2].ComovingPosition<<", v="<<Particles[2].PhysicalVelocity<<", m="<<Particles[2].Mass<<endl;
//   cout<<" Particle[end]: x="<<Particles.back().ComovingPosition<<", v="<<Particles.back().PhysicalVelocity<<", m="<<Particles.back().Mass<<endl;
}

bool IsGadget4Group(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 10)=="gadget4hdf";
}
}

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
#include <list>

#include "../snapshot.h"
#include "../mymath.h"
#include "../hdf_wrapper.h"
#include "apostle_io.h"

void create_ApostleHeader_MPI_type(MPI_Datatype& dtype)
{
  /*to create the struct data type for communication*/	
  ApostleHeader_t p;
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
void ApostleReader_t::GetParticleCountInFile(hid_t file, int np[])
{
  ReadAttribute(file, "Header", "NumPart_ThisFile", H5T_NATIVE_INT, np);
#ifdef DM_ONLY
  for(int i=0;i<TypeMax;i++)
	if(i!=TypeDM) np[i]=0;
#endif
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

void ApostleReader_t::ReadGroupParticles(int ifile, ParticleHost_t *ParticlesInFile, bool FlagReadParticleId)
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
		ParticlesThisType[i].HostId=(id[i]<0?-id[i]:id[i]);//negative means unbound 
	}
	
	H5Gclose(particle_data);
  }
  
  H5Fclose(file);
}

void ApostleReader_t::LoadSnapshot(MpiWorker_t &world, int snapshotId, vector <Particle_t> &Particles, Cosmology_t &Cosmology)
{
  SetSnapshot(snapshotId);
  
  const int root=0;
  if(world.rank()==root)
  {
    ReadHeader(0, Header);
    CompileFileOffsets(Header.NumberOfFiles);
  }
  MPI_Bcast(&Header, 1, MPI_ApostleHeader_t, root, world.Communicator);
  world.SyncContainer(np_file, MPI_HBT_INT, root);
  world.SyncContainer(offset_file, MPI_HBT_INT, root);
  
  Cosmology.Set(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);
#ifdef DM_ONLY
//   Cosmology.ParticleMass=Header.mass[TypeDM];
#endif
  
  int nfiles_skip, nfiles_end;
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

inline bool CompParticleHost(const ParticleHost_t &a, const ParticleHost_t &b)
{
  return a.HostId<b.HostId;
}

void ApostleReader_t::LoadGroups(MpiWorker_t &world, int snapshotId, vector< Halo_t >& Halos)
{//read in particle properties at the same time, to avoid particle look-up at later stage.
  SetSnapshot(snapshotId);
  
  const int root=0;
  if(world.rank()==root)
  {
    ReadHeader(0, Header);
    CompileFileOffsets(Header.NumberOfFiles);
  }
  MPI_Bcast(&Header, 1, MPI_ApostleHeader_t, root, world.Communicator);
  world.SyncContainer(np_file, MPI_HBT_INT, root);
  world.SyncContainer(offset_file, MPI_HBT_INT, root);
  
  vector <ParticleHost_t> ParticleHosts;
  int nfiles_skip, nfiles_end;
  AssignTasks(world.rank(), world.size(), Header.NumberOfFiles, nfiles_skip, nfiles_end);
  {
    HBTInt np=0;
    np=accumulate(np_file.begin()+nfiles_skip, np_file.begin()+nfiles_end, np);
    ParticleHosts.resize(np);
  }
  bool FlagReadId=true; //!HBTConfig.GroupLoadedIndex;
  
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
		ReadGroupParticles(iFile, ParticleHosts.data()+offset_file[iFile]-offset_file[nfiles_skip], FlagReadId);
	}
  }
  
  const int NullGroupId=1<<30; //1073741824
  sort(ParticleHosts.begin(), ParticleHosts.end(), CompParticleHost);
  if(!ParticleHosts.empty())
  {
    assert(ParticleHosts.back().HostId<=NullGroupId);//max haloid==NullGroupId
    assert(ParticleHosts.front().HostId>=0);//min haloid>=0
  }
  
  struct HaloLen_t
  {
    HBTInt haloid;
    HBTInt np;
    HaloLen_t()=default;
    HaloLen_t(HBTInt haloid, HBTInt np): haloid(haloid), np(np)
    {
    }
  };
  vector <HaloLen_t> HaloLen;
  
  HBTInt curr_host_id=NullGroupId;
  for(auto &&p: ParticleHosts)
  {
    if(p.HostId==NullGroupId) break;//NullGroupId comes last
    if(p.HostId!=curr_host_id)
    {
      curr_host_id=p.HostId;
      HaloLen.emplace_back(curr_host_id, 1);
    }
    else
      HaloLen.back().np++;
  }
  Halos.resize(HaloLen.size());
  for(HBTInt i=0;i<Halos.size();i++)
  {
    Halos[i].HaloId=HaloLen[i].haloid;
    Halos[i].Particles.resize(HaloLen[i].np);
  }
  auto p_in=ParticleHosts.begin();
  for(auto &&h: Halos)
  {
    for(auto &&p: h.Particles)
    {
      p=*p_in;
      ++p_in;
    }
  }
  
  VectorFree(ParticleHosts);
  
  ExchangeAndMerge(world, Halos);
  
//   cout<<Halos.size()<<" groups loaded";
//   if(Halos.size()) cout<<" : "<<Halos[0].Particles.size();
//   if(Halos.size()>1) cout<<","<<Halos[1].Particles.size()<<"...";
//   cout<<endl;
  
//   HBTInt np=0;
//   for(auto &&h: Halos)
//     np+=h.Particles.size();
//   MPI_Allreduce(MPI_IN_PLACE, &np, 1, MPI_HBT_INT, MPI_SUM, world.Communicator);
//   return np;
}

bool IsApostleGroup(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 7)=="apostle";
}

struct HaloInfo_t
{
  HBTInt id;
  HBTReal m;
  HBTxyz x;
  int order;
};
static void create_MPI_HaloInfo_t(MPI_Datatype &dtype)
{
  HaloInfo_t p;
  #define NumAttr 13
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint   offsets[NumAttr], origin,extent;
  
  MPI_Get_address(&p,&origin);
  MPI_Get_address((&p)+1,&extent);//to get the extent of s
  extent-=origin;
  
  int i=0;
  #define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
  RegisterAttr(id, MPI_HBT_INT, 1)
  RegisterAttr(m, MPI_HBT_REAL, 1)
  RegisterAttr(x[0], MPI_HBT_REAL, 3)
  RegisterAttr(order, MPI_INT, 1)
  #undef RegisterAttr
  assert(i<=NumAttr);
  
  MPI_Type_create_struct(i,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef NumAttr
}

inline bool CompHaloInfo_Id(const HaloInfo_t &a, const HaloInfo_t &b)
{
  return a.id<b.id;
}
inline bool CompHaloInfo_Order(const HaloInfo_t &a, const HaloInfo_t &b)
{
  return a.order<b.order;
}
inline bool CompHaloId(const Halo_t &a, const Halo_t &b)
{
  return a.HaloId<b.HaloId;
}
double ReduceHaloPosition(vector <HaloInfo_t>::iterator it_begin, vector <HaloInfo_t>::iterator it_end, HBTxyz &x)
{
  HBTInt i,j;
  double sx[3],origin[3],msum;
  
  if(it_begin==it_end) return 0.;
  if(it_begin+1==it_end) 
  {
    copyHBTxyz(x, it_begin->x);
    return it_begin->m;
  }
  
  sx[0]=sx[1]=sx[2]=0.;
  msum=0.;
  if(HBTConfig.PeriodicBoundaryOn)
    for(j=0;j<3;j++)
	  origin[j]=it_begin->x[j];
  
  for(auto it=it_begin;it!=it_end;++it)
  {
    HBTReal m=it->m;
    msum+=m;
    for(j=0;j<3;j++)
    if(HBTConfig.PeriodicBoundaryOn)
	    sx[j]+=NEAREST(it->x[j]-origin[j])*m;
    else
	    sx[j]+=it->x[j]*m;
  }
  
  for(j=0;j<3;j++)
  {
	  sx[j]/=msum;
	  if(HBTConfig.PeriodicBoundaryOn) 
	  {
	    sx[j]+=origin[j];
	    x[j]=position_modulus(sx[j], HBTConfig.BoxSize);
	  }
	  else
	    x[j]=sx[j];
  }
  return msum;
}
void ReduceHaloRank(vector <HaloInfo_t>::iterator it_begin, vector <HaloInfo_t>::iterator it_end, HBTxyz &step, vector <int> &dims)
{
  HBTxyz x;
  ReduceHaloPosition(it_begin, it_end,x);
  int rank=AssignCell(x, step, dims);
  for(auto it=it_begin;it!=it_end;++it)
    it->id=rank; //store destination rank in id.
}
static void DecideTargetProcessor(MpiWorker_t& world, vector< Halo_t >& Halos, vector <IdRank_t> &TargetRank)
{
  int this_rank=world.rank();
  for(auto &&h: Halos)
    h.Mass=AveragePosition(h.ComovingAveragePosition, h.Particles.data(), h.Particles.size());
   
  vector <HaloInfo_t> HaloInfoSend(Halos.size()), HaloInfoRecv;
  for(HBTInt i=0;i<Halos.size();i++)
  {
    HaloInfoSend[i].id=Halos[i].HaloId;
    HaloInfoSend[i].m=Halos[i].Mass;
    HaloInfoSend[i].x=Halos[i].ComovingAveragePosition;
//     HaloInfoSend[i].rank=this_rank;
  }
  HBTInt MaxHaloId=0; 
  if(Halos.size()) MaxHaloId=Halos.back().HaloId;
  MPI_Allreduce(MPI_IN_PLACE, &MaxHaloId, 1, MPI_HBT_INT, MPI_MAX, world.Communicator);
  HBTInt ndiv=(++MaxHaloId)/world.size();
  if(MaxHaloId%world.size()) ndiv++;
  vector <int> SendSizes(world.size(),0), SendOffsets(world.size()), RecvSizes(world.size()), RecvOffsets(world.size());
  for(HBTInt i=0;i<Halos.size();i++)
  {
    int idiv=Halos[i].HaloId/ndiv;
    SendSizes[idiv]++;
  }
  CompileOffsets(SendSizes, SendOffsets);
  MPI_Alltoall(SendSizes.data(), 1, MPI_INT, RecvSizes.data(), 1, MPI_INT, world.Communicator);
  int nhalo_recv=CompileOffsets(RecvSizes, RecvOffsets);
  HaloInfoRecv.resize(nhalo_recv);
  MPI_Datatype MPI_HaloInfo_t;
  create_MPI_HaloInfo_t(MPI_HaloInfo_t);
  MPI_Alltoallv(HaloInfoSend.data(), SendSizes.data(), SendOffsets.data(), MPI_HaloInfo_t, HaloInfoRecv.data(), RecvSizes.data(), RecvOffsets.data(), MPI_HaloInfo_t, world.Communicator);
  for(int i=0;i<nhalo_recv;i++)
    HaloInfoRecv[i].order=i;
  sort(HaloInfoRecv.begin(), HaloInfoRecv.end(), CompHaloInfo_Id);
  list <int> haloid_offsets;
  HBTInt curr_id=-1;
  for(int i=0;i<nhalo_recv;i++)
  {
    if(curr_id!=HaloInfoRecv[i].id)
    {
      haloid_offsets.push_back(i);
      curr_id=HaloInfoRecv[i].id;
    }
  }
  haloid_offsets.push_back(nhalo_recv);
  //combine coordinates and determine target  
  auto dims=ClosestFactors(world.size(), 3);
  HBTxyz step;
  for(int i=0;i<3;i++)
	step[i]=HBTConfig.BoxSize/dims[i];
  auto it_end=haloid_offsets.end(); --it_end;
  for(auto it=haloid_offsets.begin();it!=it_end;it++)
  {
    auto it_next=it;
    ++it_next;
    ReduceHaloRank(HaloInfoRecv.begin()+*it, HaloInfoRecv.begin()+*it_next, step, dims);
  }
  sort(HaloInfoRecv.begin(), HaloInfoRecv.end(), CompHaloInfo_Order);
  //send back
  MPI_Alltoallv(HaloInfoRecv.data(), RecvSizes.data(), RecvOffsets.data(), MPI_HaloInfo_t, HaloInfoSend.data(), SendSizes.data(), SendOffsets.data(), MPI_HaloInfo_t, world.Communicator);
  MPI_Type_free(&MPI_HaloInfo_t);
  
  TargetRank.resize(Halos.size());
  for(HBTInt i=0; i<TargetRank.size();i++)
  {
    TargetRank[i].Id=i;
    TargetRank[i].Rank=HaloInfoSend[i].id;
  }
  
}

void MergeHalos(vector< Halo_t >& Halos)
{
  if(Halos.empty()) return;
  sort(Halos.begin(), Halos.end(), CompHaloId);
  auto it1=Halos.begin();
  for(auto it2=it1+1;it2!=Halos.end();++it2)
  {
    if(it2->HaloId==it1->HaloId)
    {
      it1->Particles.insert(it1->Particles.end(), it2->Particles.begin(), it2->Particles.end());
    }
    else
    {
      ++it1;
      if(it2!=it1)
	*it1=move(*it2);
    }
  }
  Halos.resize(it1-Halos.begin()+1); 
  for(auto &&h: Halos)
    h.AverageCoordinates();
}

#include "../halo_particle_iterator.h"
static void ExchangeHalos(MpiWorker_t& world, vector <Halo_t>& InHalos, vector<Halo_t>& OutHalos, MPI_Datatype MPI_Halo_Shell_Type)
{
  typedef typename vector <Halo_t>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;
   
  vector <IdRank_t>TargetRank(InHalos.size());
  DecideTargetProcessor(world, InHalos, TargetRank);

  //distribute halo shells
	vector <int> SendHaloCounts(world.size(),0), RecvHaloCounts(world.size()), SendHaloDisps(world.size()), RecvHaloDisps(world.size());
	sort(TargetRank.begin(), TargetRank.end(), CompareRank);
	vector <Halo_t> InHalosSorted(InHalos.size());
	vector <HBTInt> InHaloSizes(InHalos.size());
	for(HBTInt haloid=0;haloid<InHalos.size();haloid++)
	{
	  InHalosSorted[haloid]=move(InHalos[TargetRank[haloid].Id]);
	  SendHaloCounts[TargetRank[haloid].Rank]++;
	  InHaloSizes[haloid]=InHalosSorted[haloid].Particles.size();
	}
	MPI_Alltoall(SendHaloCounts.data(), 1, MPI_INT, RecvHaloCounts.data(), 1, MPI_INT, world.Communicator);
	CompileOffsets(SendHaloCounts, SendHaloDisps);
	HBTInt NumNewHalos=CompileOffsets(RecvHaloCounts, RecvHaloDisps);
	OutHalos.resize(OutHalos.size()+NumNewHalos);
	auto NewHalos=OutHalos.end()-NumNewHalos;
	MPI_Alltoallv(InHalosSorted.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_Halo_Shell_Type, &NewHalos[0], RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_Halo_Shell_Type, world.Communicator);
  //resize receivehalos
	vector <HBTInt> OutHaloSizes(NumNewHalos);
	MPI_Alltoallv(InHaloSizes.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HBT_INT, OutHaloSizes.data(), RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HBT_INT, world.Communicator);
	for(HBTInt i=0;i<NumNewHalos;i++)
	  NewHalos[i].Particles.resize(OutHaloSizes[i]);
	
	{
	//distribute halo particles
	MPI_Datatype MPI_HBT_Particle;
	Particle_t().create_MPI_type(MPI_HBT_Particle);
	//create combined iterator for each bunch of haloes
	vector <ParticleIterator_t> InParticleIterator(world.size());
	vector <ParticleIterator_t> OutParticleIterator(world.size());
	for(int rank=0;rank<world.size();rank++)
	{
	  InParticleIterator[rank].init(InHalosSorted.begin()+SendHaloDisps[rank], InHalosSorted.begin()+SendHaloDisps[rank]+SendHaloCounts[rank]);
	  OutParticleIterator[rank].init(NewHalos+RecvHaloDisps[rank], NewHalos+RecvHaloDisps[rank]+RecvHaloCounts[rank]);
	}
	vector <HBTInt> InParticleCount(world.size(),0);
	for(HBTInt i=0;i<InHalosSorted.size();i++)
	  InParticleCount[TargetRank[i].Rank]+=InHalosSorted[i].Particles.size();
	
	MyAllToAll<Particle_t, ParticleIterator_t, ParticleIterator_t>(world, InParticleIterator, InParticleCount, OutParticleIterator, MPI_HBT_Particle);
	
	MPI_Type_free(&MPI_HBT_Particle);
	}
}

void ApostleReader_t::ExchangeAndMerge(MpiWorker_t& world, vector< Halo_t >& Halos)
{
  vector <Halo_t> LocalHalos;
  MPI_Datatype MPI_Halo_Shell_t;
  create_MPI_Halo_Id_type(MPI_Halo_Shell_t);
  ExchangeHalos(world, Halos, LocalHalos,  MPI_Halo_Shell_t);
  MPI_Type_free(&MPI_Halo_Shell_t);
  Halos.swap(LocalHalos);
  MergeHalos(Halos);
}
#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "../datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"
#include "../gravity/tree.h"
#include "../tracks.h"
	
void SubHalo_t::UpdateTrack(HBTInt snapshot_index)
{
  if(TrackId==SpecialConst::NullTrackId) return;
  
  if(0==Rank) SnapshotIndexOfLastIsolation=snapshot_index;
  if(Nbound>=LastMaxMass) 
  {
	SnapshotIndexOfLastMaxMass=snapshot_index;
	LastMaxMass=Nbound;
  }
}

HBTReal SubHalo_t::KineticDistance(const Halo_t &halo, const Snapshot_t &snapshot)
{
  HBTReal dx=PeriodicDistance(halo.ComovingPosition, ComovingPosition);
  HBTReal dv=distance(halo.PhysicalVelocity, PhysicalVelocity);
  HBTReal d=dv+snapshot.Header.Hz*snapshot.Header.ScaleFactor*dx;
  return (d>0?d:-d);
}
struct ParticleEnergy_t
{
  HBTInt pid;
  float E;
};
inline bool CompEnergy(const ParticleEnergy_t & a, const ParticleEnergy_t & b)
{
  return (a.E<b.E);
};
static HBTInt PartitionBindingEnergy(vector <ParticleEnergy_t> &Elist, const size_t len)
/*sort Elist to move unbound particles to the end*/
{//similar to the C++ partition() func
  if(len==0) return 0;
  if(len==1) return Elist[0].E<0;
  
  ParticleEnergy_t Etmp=Elist[0]; 
  auto iterforward=Elist.begin(), iterbackward=Elist.begin()+len;
  while(true)
  {
	//iterforward is a void now, can be filled
	while(true)
	{
	  iterbackward--;
	  if(iterbackward==iterforward)
	  {
		*iterforward=Etmp;
		if(Etmp.E<0) iterbackward++;
		return iterbackward-Elist.begin();
	  }
	  if(iterbackward->E<0) break;
	}
	*iterforward=*iterbackward;
	//iterbackward is a void now, can be filled
	while(true)
	{
	  iterforward++;
	  if(iterforward==iterbackward)
	  {
		*iterbackward=Etmp;
		if(Etmp.E<0) iterbackward++;
		return iterbackward-Elist.begin();
	  }
	  if(iterforward->E>0) break;
	}
	*iterbackward=*iterforward;
  }
}
static void PopMostBoundParticle(ParticleEnergy_t * Edata, const HBTInt Nbound)
{
  HBTInt imin=0;
  for(HBTInt i=1;i<Nbound;i++)
  {
	if(Edata[i].E<Edata[imin].E) imin=i;
  }
  if(imin!=0) swap(Edata[imin], Edata[0]);
}
void SubHalo_t::Unbind(const Snapshot_t &snapshot)
{//the reference frame should already be initialized before unbinding.
  if(1==Particles.size()) return;
  
  HBTInt OldMostBoundParticle=Particles[0];
  OctTree_t tree;
  tree.Reserve(Particles.size());
  Nbound=Particles.size(); //start from full set
  HBTInt Nlast; 
  
  vector <ParticleEnergy_t> Elist(Nbound);
	for(HBTInt i=0;i<Nbound;i++)
	  Elist[i].pid=Particles[i];
	while(true)
	{
		Nlast=Nbound;
		tree.Build(Nlast, Particles.data(), snapshot);
	  #pragma omp parallel for if(Nlast>100)
	  for(HBTInt i=0;i<Nlast;i++)
	  {
		HBTInt pid=Elist[i].pid;
		Elist[i].E=tree.BindingEnergy(snapshot.GetComovingPosition(pid), snapshot.GetPhysicalVelocity(pid), ComovingPosition, PhysicalVelocity, snapshot.GetParticleMass(pid));
	  }
		Nbound=PartitionBindingEnergy(Elist, Nlast);//TODO: parallelize this.
		if(Nbound<HBTConfig.MinNumPartOfSub)
		{
		  Nbound=0;
		  Elist[0].pid=OldMostBoundParticle;
		}
		else
		{
		  sort(Elist.begin()+Nbound, Elist.begin()+Nlast, CompEnergy); //only sort the unbound part
		  PopMostBoundParticle(Elist.data(), Nbound);
		}
	  if(Nbound<HBTConfig.MinNumPartOfSub)
	  {
		for(HBTInt i=0;i<Nlast;i++) Particles[i]=Elist[i].pid;//copy the sorted particles this round. ToDo: optimize this away (passing Elist to tree? maybe pass functions to tree, or a ParticleList/minisnapshot base class with GetPos,GetMass etc.)
	  }
	  {
		copyHBTxyz(ComovingPosition, snapshot.GetComovingPosition(Elist[0].pid));
		copyHBTxyz(PhysicalVelocity, snapshot.GetPhysicalVelocity(Elist[0].pid));
	  }
		if(0==Nbound||Nbound>Nlast*HBTConfig.BoundMassPrecision)  break;
	}
	if(Nbound)
	{
	  sort(Elist.begin(), Elist.begin()+Nbound, CompEnergy); //sort the self-bound part
	  for(HBTInt i=0;i<Nbound;i++) Particles[i]=Elist[i].pid;//unbound part already copied inside the loop.
		Nlast=Particles.size();
		if(Nlast>Nbound*HBTConfig.SourceSubRelaxFactor) Nlast=Nbound*HBTConfig.SourceSubRelaxFactor;
	}
	else
	{
	  Nbound=1;
	  Nlast=1;//what if this is a central?? any fixes?
	}
  Particles.resize(Nlast);
  //todo: output angular momentum and total energy as well, for calculation of spin.
}
void MemberShipTable_t::ResizeAllMembers(size_t n)
{
  Mem_AllMembers.resize(n);
  size_t offset=Mem_AllMembers.data()-AllMembers.data();
  if(offset)
  {
	AllMembers.Bind(n, AllMembers.data()+offset);
	for(HBTInt i=0;i<Mem_SubGroups.size();i++)
	  Mem_SubGroups[i].Bind(Mem_SubGroups[i].data()+offset);
  }
}
void MemberShipTable_t::Init(const HBTInt nhalos, const HBTInt nsubhalos, const float alloc_factor)
{
  Mem_SubGroups.clear();
  Mem_SubGroups.resize(nhalos+1);
  SubGroups.Bind(nhalos, Mem_SubGroups.data()+1);

  Mem_AllMembers.clear();
  Mem_AllMembers.reserve(nsubhalos*alloc_factor); //allocate more for seed haloes.
  Mem_AllMembers.resize(nsubhalos);
  AllMembers.Bind(nsubhalos, Mem_AllMembers.data());
}
void MemberShipTable_t::BindMemberLists()
{
  HBTInt offset=0;
  for(HBTInt i=0;i<Mem_SubGroups.size();i++)
  {
	Mem_SubGroups[i].Bind(Mem_SubGroups[i].size(), &(Mem_AllMembers[offset]));
	offset+=Mem_SubGroups[i].size();
	Mem_SubGroups[i].ReBind(0);
  }
}
void MemberShipTable_t::CountMembers(const SubHaloList_t& SubHalos)
{//todo: parallelize this..
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
	SubGroups[SubHalos[subid].HostHaloId].IncrementBind();
}
void MemberShipTable_t::FillMemberLists(const SubHaloList_t& SubHalos)
{
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
  {
	SubGroups[SubHalos[subid].HostHaloId].PushBack(subid);
  }
}
struct CompareMass_t
{
  const SubHaloList_t * SubHalos;
  CompareMass_t(const SubHaloList_t & subhalos)
  {
	SubHalos=&subhalos;
  }
  bool operator () (const HBTInt & i, const HBTInt & j)
  {
	return (*SubHalos)[i].Nbound>(*SubHalos)[j].Nbound;
  }
};
void MemberShipTable_t::SortMemberLists(const SubHaloList_t & SubHalos)
{
  CompareMass_t compare_mass(SubHalos);
#pragma omp for
  for(HBTInt i=-1;i<SubGroups.size();i++)
	std::sort(SubGroups[i].begin(), SubGroups[i].end(), compare_mass);
}
void MemberShipTable_t::SortSatellites(const SubHaloList_t & SubHalos)
/*central subhalo not changed*/
{
  CompareMass_t compare_mass(SubHalos);
  for(HBTInt i=0;i<SubGroups.size();i++)
	std::sort(SubGroups[i].begin()+1, SubGroups[i].end(), compare_mass);
}
void MemberShipTable_t::AssignRanks(SubHaloList_t& SubHalos)
{
#pragma omp for
  for(HBTInt haloid=-1;haloid<SubGroups.size();haloid++)
  {
	MemberList_t & SubGroup=SubGroups[haloid];
	for(HBTInt i=0;i<SubGroup.size();i++)
	  SubHalos[SubGroup[i]].Rank=i;
  }
}
void MemberShipTable_t::CountBirth()
{
static HBTInt nbirth;
#pragma omp single
nbirth=0;
#pragma omp for reduction(+:nbirth)
  for(HBTInt hostid=0;hostid<SubGroups.size();hostid++)
	if(SubGroups[hostid].size()==0)  nbirth++;
#pragma omp single
NBirth=nbirth;	
}
/*
inline bool SubHaloSnapshot_t::CompareHostAndMass(const HBTInt& subid_a, const HBTInt& subid_b)
{//ascending in host id, descending in mass inside each host, and put NullHaloId to the beginning.
  SubHalo_t a=SubHalos[subid_a], b=SubHalos[subid_b];
  
  if(a.HostHaloId==b.HostHaloId) return (a.Nbound>b.Nbound);
  
  return (a.HostHaloId<b.HostHaloId); //(a.HostHaloId!=SpecialConst::NullHaloId)&&
}*/
void MemberShipTable_t::Build(const HBTInt nhalos, const SubHaloList_t & SubHalos)
{
  #pragma omp single
  {
  Init(nhalos, SubHalos.size());
  CountMembers(SubHalos);
  BindMemberLists();
  FillMemberLists(SubHalos);
  }
  SortMemberLists(SubHalos);
  CountBirth();
//   std::sort(AllMembers.begin(), AllMembers.end(), CompareHostAndMass);
}
void MemberShipTable_t::SubIdToTrackId(const SubHaloList_t& SubHalos)
{
  /*not necessary currently*/
   for(HBTInt i=0;i<AllMembers.size();i++)
	AllMembers[i]=SubHalos[AllMembers[i]].TrackId;
}
void MemberShipTable_t::TrackIdToSubId(SubHaloList_t& SubHalos)
{
cout<<"Warning: TrackIdToSubId ToBe fully Implemented!\n";
// exit(1);
}

void SubHaloSnapshot_t::ParticleIdToIndex(const Snapshot_t& snapshot)
{//also bind to snapshot
#pragma omp single
  SnapshotPointer=&snapshot;
#pragma omp for
	for(HBTInt subid=0;subid<SubHalos.size();subid++)
	{
	  SubHalo_t::ParticleList_t & Particles=SubHalos[subid].Particles;
	  HBTInt nP=Particles.size();
	  for(HBTInt pid=0;pid<nP;pid++)
		Particles[pid]=snapshot.GetParticleIndex(Particles[pid]);
	}
}
void SubHaloSnapshot_t::ParticleIndexToId()
{
#pragma omp for
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
  {
	SubHalo_t::ParticleList_t & Particles=SubHalos[subid].Particles;
	HBTInt nP=Particles.size();
	for(HBTInt pid=0;pid<nP;pid++)
	  Particles[pid]=SnapshotPointer->GetParticleId(Particles[pid]);
  }
#pragma omp single
  SnapshotPointer=nullptr;
}
void SubHaloSnapshot_t::BuildHDFDataType()
{
  hsize_t dims=3;
  H5::ArrayType H5T_HBTxyz(H5T_HBTReal, 1, &dims);
  #define InsertMember(x,t) H5T_SubHaloInMem.insertMember(#x, HOFFSET(SubHalo_t, x), t)//;cout<<#x<<": "<<HOFFSET(SubHalo_t, x)<<endl
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
  H5T_SubHaloInDisk.copy(H5T_SubHaloInMem);
  H5T_SubHaloInDisk.pack(); //clear fields not added.
//   SubHalo_t s;
//   cout<<(char *)&s.TrackId-(char *)&s<<","<<(char *)&s.Nbound-(char *)&s<<","<<(char *)&s.ComovingPosition-(char *)&s<<","<<(char *)&s.Particles-(char *)&s<<endl;
}
void SubHaloSnapshot_t::GetSubFileName(string &filename)
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SubSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
  filename=formater.str();
}
void SubHaloSnapshot_t::GetSrcFileName(string &filename)
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SrcSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
  filename=formater.str();
}
void SubHaloSnapshot_t::Load(int snapshot_index, bool load_src)
{
  if(snapshot_index<HBTConfig.MinSnapshotIndex)
  {
	cout<<"Warning: snapshot index "<<snapshot_index<<" is below MinSnapshotIndex of "<<HBTConfig.MinSnapshotIndex;
	cout<<", I wil not load anything.\n";
	return;
  }
  
  SetSnapshotIndex(snapshot_index);
  string filename;
  GetSubFileName(filename);
  H5::H5File file(filename.c_str(), H5F_ACC_RDONLY);
  
  H5::DataSet dset(file.openDataSet("SnapshotId"));
  HBTInt snapshot_id;
  dset.read(&snapshot_id, H5T_HBTInt);//HBTInt does not have to be the same as the datatype in the file.
  dset.close();
  assert(snapshot_id==SnapshotId);
  
  hsize_t dims[1];
  dset=file.openDataSet("SubHalos");
  dset.getSpace().getSimpleExtentDims(dims);
  HBTInt nsubhalos=dims[0];
  SubHalos.resize(nsubhalos);
  if(nsubhalos)	dset.read(SubHalos.data(), H5T_SubHaloInMem);
  dset.close();
  
  file.openDataSet("/Membership/NumberOfNewSubHalos").read(&MemberTable.NBirth, H5T_HBTInt);
  file.openDataSet("/Membership/NumberOfNewSubHalos").read(&MemberTable.NFake, H5T_HBTInt);
  dset=file.openDataSet("/Membership/GroupedTrackIds");
  dset.getSpace().getSimpleExtentDims(dims);
  HBTInt nhalos=dims[0]-1;
  vector <hvl_t> vl(dims[0]);
  H5::VarLenType H5T_HBTIntArr(&H5T_HBTInt);
  if(nsubhalos)
  {
	dset.read(vl.data(), H5T_HBTIntArr);
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
	dset.vlenReclaim(vl.data(), H5T_HBTIntArr, dset.getSpace());
	assert(offset==nsubhalos);
  }
  dset.close();
  
  if(0==nsubhalos) return;
  
  vl.resize(nsubhalos);
  if(!load_src)
  {
	dset=file.openDataSet("SubHaloParticles");
	dset.getSpace().getSimpleExtentDims(dims);
	assert(dims[0]==nsubhalos);
	dset.read(vl.data(), H5T_HBTIntArr);
	for(HBTInt i=0;i<nsubhalos;i++)
	{
	  SubHalos[i].Particles.resize(vl[i].len);
	  memcpy(SubHalos[i].Particles.data(), vl[i].p, sizeof(HBTInt)*vl[i].len);
	}
	dset.vlenReclaim(vl.data(), H5T_HBTIntArr,dset.getSpace());
	dset.close();
  }
  else
  {
	GetSrcFileName(filename);
	H5::H5File file2(filename.c_str(), H5F_ACC_RDONLY);
	dset=file2.openDataSet("SrcHaloParticles");
	dset.getSpace().getSimpleExtentDims(dims);
	assert(dims[0]==nsubhalos);
	dset.read(vl.data(), H5T_HBTIntArr);
	for(HBTInt i=0;i<nsubhalos;i++)
	{
	  SubHalos[i].Particles.resize(vl[i].len);
	  memcpy(SubHalos[i].Particles.data(), vl[i].p, sizeof(HBTInt)*vl[i].len);
	}
	dset.vlenReclaim(vl.data(), H5T_HBTIntArr,dset.getSpace());
	dset.close();
  }
  cout<<SubHalos.size()<<" subhaloes loaded at snapshot "<<SnapshotIndex<<"("<<SnapshotId<<")\n";
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
void SubHaloSnapshot_t::Save()
{
  stringstream formater;
  formater<<HBTConfig.SubhaloPath<<"/SubSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
  string filename(formater.str());
  cout<<"Saving to "<<filename<<"..."<<endl;
  H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);

  hsize_t ndim=1, dim_atom[]={1};
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5::PredType::NATIVE_INT);
  
  hsize_t dim_sub[]={SubHalos.size()};
  writeHDFmatrix(file, SubHalos.data(), "SubHalos", ndim, dim_sub, H5T_SubHaloInMem, H5T_SubHaloInDisk); 
  
  H5::Group datagrp(file.createGroup("/Membership"));
  H5LTset_attribute_string(file.getId(),"/Membership","Comment","List of subhaloes in each group.");
  writeHDFmatrix(datagrp, &MemberTable.NBirth, "NumberOfNewSubHalos", ndim, dim_atom, H5T_HBTInt);
  writeHDFmatrix(datagrp, &MemberTable.NFake, "NumberOfFakeHalos", ndim, dim_atom, H5T_HBTInt);
 
  MemberTable.SubIdToTrackId(SubHalos);
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
  writeHDFmatrix(datagrp, vl.data(), "GroupedTrackIds", ndim, dim_grp, H5T_HBTIntArr);
  H5LTset_attribute_string(datagrp.getId(),"GroupedTrackIds","Comment","Nhalo+1 groups. The last group contain tracks outside any host halo (i.e., field subhaloes).");
  
  //now write the particle list for each subhalo
  vl.resize(SubHalos.size());
  for(HBTInt i=0;i<vl.size();i++)
  {
	vl[i].len=SubHalos[i].Nbound;
	vl[i].p=SubHalos[i].Particles.data();
  }
  writeHDFmatrix(file, vl.data(), "SubHaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  
  file.close();
  
  formater.str("");
  formater.clear();
  formater<<HBTConfig.SubhaloPath<<"/SrcSnap_"<<setw(3)<<setfill('0')<<SnapshotIndex<<".hdf5"; //or use snapshotid
  filename=formater.str();
  cout<<"Saving to "<<filename<<"..."<<endl;
  file=H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
  writeHDFmatrix(file, &SnapshotId, "SnapshotId", ndim, dim_atom, H5::PredType::NATIVE_INT);
  for(HBTInt i=0;i<vl.size();i++)
	vl[i].len=SubHalos[i].Particles.size();
  writeHDFmatrix(file, vl.data(), "SrcHaloParticles", ndim, dim_sub, H5T_HBTIntArr);
  file.close();
  cout<<SubHalos.size()<<" subhaloes saved: "<<MemberTable.NBirth<<" birth, "<< MemberTable.NFake<<" fake.\n";
}
void SubHaloSnapshot_t::AssignHosts(const HaloSnapshot_t &halo_snap)
{
  static vector<HBTInt> ParticleToHost;//to make it shared
#pragma omp single
{
   ParticleToHost.assign(SnapshotPointer->GetNumberOfParticles(), SpecialConst::NullHaloId);
   ParallelizeHaloes=halo_snap.NumPartOfLargestHalo<0.1*halo_snap.TotNumberOfParticles;//no dominating objects
}
#pragma omp for
  for(HBTInt haloid=0;haloid<halo_snap.Halos.size();haloid++)
  {
	const Halo_t::ParticleList_t & Particles=halo_snap.Halos[haloid].Particles;
	for(HBTInt i=0;i<Particles.size();i++)
	  ParticleToHost[Particles[i]]=haloid;
  }
#pragma omp for 
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
  {
	//rely on most-bound particle
	SubHalos[subid].HostHaloId=ParticleToHost[SubHalos[subid].Particles[0]];
	//alternatives: CoreTrack; Split;
  }
#pragma omp single
  vector <HBTInt>().swap(ParticleToHost);//free the memory.
#pragma omp single nowait
  if(HBTConfig.TrimNonHostParticles)
  {
	cout<<"Error: TrimNonHostParticles not implemented yet...\n";
	exit(1);
  }
  MemberTable.Build(halo_snap.Halos.size(), SubHalos);
//   MemberTable.AssignRanks(SubHalos); //not needed here
}
inline HBTInt GetCoreSize(HBTInt nbound)
{
  int coresize=nbound*HBTConfig.SubCoreSizeFactor;
  if(coresize<HBTConfig.SubCoreSizeMin) coresize=HBTConfig.SubCoreSizeMin;
  if(coresize>nbound) coresize=nbound;
  return coresize;
}
void SubHaloSnapshot_t::AverageCoordinates()
{
#pragma omp for
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
  {
	int coresize=GetCoreSize(SubHalos[subid].Nbound);
	SnapshotPointer->AveragePosition(SubHalos[subid].ComovingPosition, SubHalos[subid].Particles.data(), coresize);
	SnapshotPointer->AverageVelocity(SubHalos[subid].PhysicalVelocity, SubHalos[subid].Particles.data(), coresize);
  }
}

void SubHaloSnapshot_t::DecideCentrals(const HaloSnapshot_t &halo_snap)
/* to select central subhalo according to KineticDistance, and move each central to the beginning of each list in MemberTable*/
{
#pragma omp for
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
  {
	MemberShipTable_t::MemberList_t &List=MemberTable.SubGroups[hostid];
	if(List.size()>1)
	{
	  int n_major;
	  HBTInt MassLimit=SubHalos[List[0]].Nbound*HBTConfig.MajorProgenitorMassRatio;
	  for(n_major=1;n_major<List.size();n_major++)
		if(SubHalos[List[n_major]].Nbound<MassLimit) break;
	  if(n_major>1)
	  {
		HBTReal dmin=SubHalos[List[0]].KineticDistance(halo_snap.Halos[hostid], *SnapshotPointer);
		int icenter=0;
		for(int i=1;i<n_major;i++)
		{
		  HBTReal d=SubHalos[List[i]].KineticDistance(halo_snap.Halos[hostid], *SnapshotPointer);
		  if(dmin>d)
		  {
			dmin=d;
			icenter=i;
		  }
		}
		if(icenter)  swap(List[0], List[icenter]);
	  }
	}
  }
}

void SubHaloSnapshot_t::FeedCentrals(HaloSnapshot_t& halo_snap)
/* replace centrals with host particles;
 * create a new central if there is none
 * initialize new central with host halo center coordinates
 * halo_snap is rendered unspecified upon return (its particles have been swapped to subhaloes).
 */
{
  static HBTInt Npro;
  #pragma omp single
  {
  Npro=SubHalos.size();
  //SubHalos.reserve(Snapshot->GetNumberOfParticles()*0.1);//reserve enough	branches.......
  SubHalos.resize(Npro+MemberTable.NBirth);
  }
  #pragma omp for
  for(HBTInt hostid=0;hostid<halo_snap.Halos.size();hostid++)
  { 
	MemberShipTable_t::MemberList_t & Members=MemberTable.SubGroups[hostid];
	if(0==Members.size()) //create a new sub
	{
	  HBTInt subid;
	  #pragma omp critical(AddNewSub)
	  {
		subid=Npro++;
	  }
	  SubHalos[subid].HostHaloId=hostid;
	  copyHBTxyz(SubHalos[subid].ComovingPosition, halo_snap.Halos[hostid].ComovingPosition); 
	  copyHBTxyz(SubHalos[subid].PhysicalVelocity, halo_snap.Halos[hostid].PhysicalVelocity);
	  SubHalos[subid].Particles.swap(halo_snap.Halos[hostid].Particles);
	}
	else
	  SubHalos[Members[0]].Particles.swap(halo_snap.Halos[hostid].Particles); //reuse the halo particles
  }
  #pragma omp single
  halo_snap.Clear();//to avoid misuse
}
void SubHaloSnapshot_t::PrepareCentrals(HaloSnapshot_t &halo_snap)
{
  halo_snap.AverageCoordinates();
  AverageCoordinates();
  DecideCentrals(halo_snap);
  FeedCentrals(halo_snap);
}

void SubHaloSnapshot_t::RefineParticles()
{//it's more expensive to build an exclusive list. so do inclusive here. 
  //TODO: ensure the inclusive unbinding is stable (contaminating particles from big subhaloes may hurdle the unbinding
#ifdef _OPENMP
 if(ParallelizeHaloes) cout<<"Unbinding with HaloPara...\n";
 else cout<<"Unbinding with ParticlePara...\n";
#else
 cout<<"Unbinding..."<<endl;
#endif  
#pragma omp parallel for if(ParallelizeHaloes)
  for(HBTInt subid=0;subid<SubHalos.size();subid++)
  {
	try
	{
	  SubHalos[subid].Unbind(*SnapshotPointer);
	}
	catch(OctTreeExceeded_t &tree_exception)
	{
	  cerr<<"Error: "<<tree_exception.what()<<" in subhalo "<<subid<<endl;
	  exit(1);
	}
  }
}

void SubHaloSnapshot_t::RegisterNewTracks()
/*assign trackId to new bound ones, remove unbound ones, and record membership*/
{
  HBTInt NTot=SubHalos.size();
  HBTInt TrackId=NTot-MemberTable.NBirth, NFake=0;
  MemberTable.ResizeAllMembers(NTot);
  for(HBTInt i=TrackId;i<NTot;i++)
  {
	if(SubHalos[i].Nbound>1)
	{
	  if(i!=TrackId)
		SubHalos[i].MoveTo(SubHalos[TrackId]);
	  SubHalos[TrackId].TrackId=TrackId;
	  MemberTable.AllMembers[TrackId]=TrackId; //this trackId is also the subhalo index
// 	  assert(MemberTable.SubGroups[SubHalos[TrackId].HostHaloId].size()==0);
	  MemberTable.SubGroups[SubHalos[TrackId].HostHaloId].Bind(1, &MemberTable.AllMembers[TrackId]);
	  TrackId++;
	}
  }
  MemberTable.NFake=NTot-TrackId;
  MemberTable.NBirth-=MemberTable.NFake;
  SubHalos.resize(TrackId);
}
void SubHaloSnapshot_t::UpdateTracks()
{
  /*renew ranks after unbinding*/
#pragma omp single
  RegisterNewTracks();
  MemberTable.SortMemberLists(SubHalos);//reorder
  MemberTable.AssignRanks(SubHalos);
#pragma omp for
  for(HBTInt i=0;i<SubHalos.size();i++)
	SubHalos[i].UpdateTrack(SnapshotIndex);
}

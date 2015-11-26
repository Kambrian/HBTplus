#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <glob.h>
#include <climits>
#include <algorithm>
#include <chrono>

#include "../mymath.h"
#include "../halo.h"

#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
#define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)

struct GroupV4Header_t
{
  int Ngroups;
  int Nsubgroups;
  int Nids;
  int TotNgroups;
  int TotNsubgroups;
  int TotNids;
  int num_files;//long? no, but padding may exist.
  double time;
  double redshift;
  double HubbleParam;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  int flag_doubleprecision;//long? no, but padding may exist
};
static bool GetGroupFileByteOrder(const char *filename, const int FileCounts, const int GroupFileVariant)
{
  /* to check whether byteswap is needed, return true if yes, false if no, exit if error*/
  int Nfiles,n,ns;
  long offset;
  FILE *fp;
  
  if(GROUP_FORMAT_GADGET4==GroupFileVariant)
	n=sizeof(GroupV4Header_t);
  else
	n=FileCounts;
  
  ns=n;
  swap_Nbyte(&ns,1,sizeof(ns));	
  
  myfopen(fp,filename,"r");
  switch(GroupFileVariant)
  {
	case GROUP_FORMAT_GADGET4:
	  offset=0;
	  break;
	case GROUP_FORMAT_GADGET3_INT:
	case GROUP_FORMAT_GADGET3_LONG:
	  offset=3*sizeof(int)+sizeof(long long);  
	  break;
	case GROUP_FORMAT_GADGET2_INT:
	case GROUP_FORMAT_GADGET2_LONG:
	default:
	  offset=3*sizeof(int); 
  }  
  fseek(fp,offset,SEEK_SET);
  size_t tmp_size=fread(&Nfiles,sizeof(int),1,fp);
  fclose(fp);
  
  if(Nfiles==n) return false;
  if(Nfiles==ns) return true;
  
  cerr<<"endianness check failed for: "<<filename<<", file format not expected:"<<Nfiles<<';'<<n<<';'<<ns<<endl;
  exit(1);
}
int HaloSnapshot_t::CountFiles()
{
  string format;
  bool IsSubFile, NeedByteSwap;
  int FileCounts;
  GetFileNameFormat(format, FileCounts, IsSubFile, NeedByteSwap);
  return FileCounts;
}
void HaloSnapshot_t::GetFileNameFormat(string &format, int &FileCounts, bool &IsSubFile, bool &NeedByteSwap)
{
  IsSubFile=false;
  FileCounts=1;
  char filename[1024], pattern[1024], basefmt[1024], fmt[1024];
  const int ifile=0;
  sprintf(basefmt, "%s/groups_%03d/subhalo_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId,SnapshotId);
  sprintf(fmt, "%s.%%d", basefmt);
  sprintf(filename, fmt, "tab", ifile);
  if(file_exist(filename))
  {
	IsSubFile=true;
	sprintf(pattern, basefmt, "tab");
	strcat(pattern, ".*");
	FileCounts=count_pattern_files(pattern);
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts, HBTConfig.GroupFileVariant);
	format=fmt;
	return;
  }
  
  sprintf(basefmt, "%s/groups_%03d/group_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId,SnapshotId);
  sprintf(fmt, "%s.%%d", basefmt);
  sprintf(filename, fmt, "tab", ifile);
  if(file_exist(filename))
  {
	sprintf(pattern, basefmt, "tab");
	strcat(pattern, ".*");
	FileCounts=count_pattern_files(pattern);
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts, HBTConfig.GroupFileVariant);
	format=fmt;
	return;
  }
  
  sprintf(fmt, "%s/subhalo_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId);
  sprintf(filename, fmt, "tab");
  if(file_exist(filename))
  {
	IsSubFile=true;
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts, HBTConfig.GroupFileVariant);
	format=fmt;
	return;
  }
  
  sprintf(fmt, "%s/group_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId);
  sprintf(filename, fmt, "tab");
  if(file_exist(filename))
  {
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts, HBTConfig.GroupFileVariant);
	format=fmt;
	return;
  }
  
  // 	return 0;
  cerr<<"Error: no group files found under "<<HBTConfig.HaloPath<<" at snapshot "<<SnapshotId<<endl;
  exit(1);
}
//for locally loaded file contents
vector <HBTInt> HaloParticleBuffer;
vector <int> HaloLenBuffer,HaloOffsetBuffer;
#define READ_META_DATA 0
#define READ_LEN_OFFSET 1
#define READ_PARTICLES 2
struct FileAssignment_t
{
  int nfile_begin, nfile_end;
  HBTInt npart_begin, npart_end;
  HBTInt fileoffset_begin, fileoffset_end;
  HBTInt haloid_begin, nhalo;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
	ar & nfile_begin;
	ar & nfile_end;
	ar & npart_begin;
	ar & npart_end;
	ar & fileoffset_begin;
	ar & fileoffset_end;
	ar & haloid_begin;
	ar & nhalo;
  }
};
BOOST_IS_MPI_DATATYPE(FileAssignment_t)

void AssignHaloTasks(int nworkers, int npart_tot, const vector <HBTInt> &FileOffsets, int & npart_begin, FileAssignment_t &task)
{
  int npart_this=(npart_tot-npart_begin)/nworkers;
  int npart_end=npart_begin+npart_this;
  task.haloid_begin=lower_bound(HaloOffsetBuffer.begin(), HaloOffsetBuffer.end()-1, npart_begin)-HaloOffsetBuffer.begin();
  HBTInt endhalo=lower_bound(HaloOffsetBuffer.begin(), HaloOffsetBuffer.end()-1, npart_end)-HaloOffsetBuffer.begin();
  task.nhalo=endhalo-task.haloid_begin;
  task.npart_begin=npart_begin;//to subtract fileoffset!!
  task.npart_end=HaloOffsetBuffer.at(endhalo);
  
  task.nfile_begin=upper_bound(FileOffsets.begin(), FileOffsets.end(), task.npart_begin)-1-FileOffsets.begin();
  task.nfile_end=lower_bound(FileOffsets.begin(), FileOffsets.end(), task.npart_end)-FileOffsets.begin();
  
  task.fileoffset_begin=FileOffsets[task.nfile_begin];
  task.fileoffset_end=FileOffsets[task.nfile_end-1];
  
  npart_begin=task.npart_end;
}
void HaloSnapshot_t::Load(mpi::communicator & world, int snapshot_index)
{
  SetSnapshotIndex(snapshot_index);
  int FileCounts=CountFiles();

  /*distribute tasks*/
  vector <FileAssignment_t> alltasks;
  FileAssignment_t thistask;
  HBTInt Ngroups=0,Nparticles=0;
  GroupFileSize_t filesize;
  if(world.rank()==0)
  {
	vector <HBTInt> FileOffsets;
	FileOffsets.resize(FileCounts);
	for(int iFile=0;iFile<FileCounts;iFile++)
	{
	  FileOffsets[iFile]=Nparticles;
	  ReadFile(iFile, filesize, READ_META_DATA);
	  Ngroups+=filesize.NumberOfGroups;
	  Nparticles+=filesize.NumberOfParticles;
	}
	HaloLenBuffer.reserve(Ngroups);
	HaloOffsetBuffer.reserve(Ngroups);
	for(int iFile=0;iFile<FileCounts;iFile++)
	  ReadFile(iFile, filesize, READ_LEN_OFFSET);
	HaloOffsetBuffer.push_back(Nparticles);//end offset
	
	alltasks.resize(world.size());
	{
	  int nworkers=world.size(), npart_begin=0; 
	  for(int rank=0;rank<world.size();rank++)
		AssignHaloTasks(nworkers--, Nparticles, FileOffsets, npart_begin, alltasks[rank]);
	}
  }
  scatter(world, alltasks, thistask, 0);
  if(world.rank()==0)
  {
	for(int rank=1;rank<world.size();rank++)
	  world.send(rank, 0, HaloLenBuffer.data()+alltasks[rank].haloid_begin, alltasks[rank].nhalo);
	HaloLenBuffer.resize(thistask.nhalo);
  }
  else
  {
	HaloLenBuffer.resize(thistask.nhalo);
	world.recv(0,0,HaloLenBuffer.data(), HaloLenBuffer.size());
  }
  
  /* read particles*/
  HaloParticleBuffer.reserve(thistask.npart_end-thistask.npart_begin);
  for(int iFile=thistask.nfile_begin;iFile<thistask.nfile_end;iFile++)
  {
	HBTInt start_particle=0, end_particle=-1;
	if(iFile==thistask.nfile_begin) start_particle=thistask.npart_begin-thistask.fileoffset_begin;
	if(iFile==thistask.nfile_end-1) end_particle=thistask.npart_end-thistask.fileoffset_end;
	ReadFile(iFile, filesize, READ_PARTICLES, start_particle, end_particle);
  }
  
  /* populate haloes*/
  Halos.resize(thistask.nhalo);
  auto p=HaloParticleBuffer.begin();
  for(HBTInt i=0;i<Halos.size();i++)
  {
	Halos[i].Particles.assign(p, p+HaloLenBuffer[i]);
	p+=HaloLenBuffer[i];
  }
  assert(p==HaloLenBuffer.end());
  
  for(auto &&np: HaloLenBuffer)
  {
	TotNumberOfParticles+=np;//local
	if(NumPartOfLargestHalo<np) NumPartOfLargestHalo=np;//local
  }
  assert(TotNumberOfParticles==(thistask.npart_end-thistask.npart_begin));
  
  /*clear up buffers*/
  vector <HBTInt>().swap(HaloParticleBuffer);
  vector <int>().swap(HaloLenBuffer);
  vector <int>().swap(HaloOffsetBuffer);
  
  /*distribute groups into domains*/
  
  cout<<"Finished reading "<<Halos.size()<<" groups ("<<TotNumberOfParticles<<" particles) from file "<<thistask.nfile_begin<<" to "<<thistask.nfile_end-1<<" (total "<<FileCounts<<" files) on thread "<<world.rank()<<endl;  
}

void HaloSnapshot_t::ReadFile(int iFile, GroupFileSize_t &filesize, int read_level, HBTInt start_particle, HBTInt end_particle)
{
  switch(HBTConfig.GroupFileVariant)
  {
	case GROUP_FORMAT_GADGET2_INT:
	case GROUP_FORMAT_GADGET3_INT:
	  ReadFileV2V3<int>(iFile, filesize, read_level, start_particle, end_particle);
	  break;
	case GROUP_FORMAT_GADGET2_LONG:
	case GROUP_FORMAT_GADGET3_LONG:
	  ReadFileV2V3<long>(iFile, filesize, read_level, start_particle, end_particle);
	  break;
	default:
	  cerr<<"GroupFileVariant="<<HBTConfig.GroupFileVariant<<" not implemented yet.\n";
	  exit(1);
  }
}
void HaloSnapshot_t::Clear()
/* call this to reset the HaloSnapshot to empty.
 This is usually not necessary because the destructor will release the memory automatically*/
{
//   Halos.clear(); //this does not actually clear
  HaloList_t().swap(Halos);//this deeply cleans it
}

#define IS_GROUP_V3 (GROUP_FORMAT_GADGET3_INT==HBTConfig.GroupFileVariant||GROUP_FORMAT_GADGET3_LONG==HBTConfig.GroupFileVariant)

template <class PIDtype_t>
void HaloSnapshot_t::ReadFileV2V3(int iFile, GroupFileSize_t &filesize, int read_level, HBTInt start_particle, HBTInt end_particle)
/*read_level: 0: meta data only ; 1: further read group len and offset; 2: further read particle lists*/
{
  FILE *fd;
  char filename[1024];
  int Ngroups, TotNgroups, Nids, NFiles;
  long long TotNids;
  bool NeedByteSwap;
  bool IsGroupV3=IS_GROUP_V3;
  
  string filename_format;
  bool IsSubFile;
  int FileCounts;
  GetFileNameFormat(filename_format, FileCounts, IsSubFile, NeedByteSwap);
  
	if(FileCounts>1)
	  sprintf(filename, filename_format.c_str(), "tab", iFile);
	else
	  sprintf(filename, filename_format.c_str(), "tab");
		
	myfopen(fd,filename,"r");
	myfread(&Ngroups, sizeof(Ngroups), 1, fd);
	if(IsGroupV3)
	{
	  myfread(&TotNgroups, sizeof(TotNgroups), 1, fd);
	  myfread(&Nids, sizeof(Nids), 1, fd);
	  myfread(&TotNids,sizeof(TotNids),1,fd);
	}
	else
	{
	  myfread(&Nids, sizeof(Nids), 1, fd);
	  myfread(&TotNgroups, sizeof(TotNgroups), 1, fd);
	}
	filesize.NumberOfGroups=Ngroups;
	filesize.NumberOfParticles=Nids;
	myfread(&NFiles, sizeof(NFiles), 1, fd);
	if(IsSubFile)
	{
	  int Nsub,TotNsub;
	  myfread(&Nsub,sizeof(int),1,fd);
	  myfread(&TotNsub,sizeof(int),1,fd);
	}
	if(FileCounts!=NFiles)
	{
	  cout<<"File count mismatch for file "<<filename<<": expect "<<FileCounts<<", got "<<NFiles<<endl;
	  exit(1);
	}
	if(read_level==0) 
	{
	  fclose(fd);
	  return;
	}
	
	if(read_level==1)
	{
	  HaloLenBuffer.resize(HaloLenBuffer.size()+Ngroups);
	  HaloOffsetBuffer.resize(HaloOffsetBuffer.size()+Ngroups);
	  myfread(HaloLenBuffer.data()+HaloLenBuffer.size()-Ngroups, sizeof(int), Ngroups, fd);
	  myfread(HaloOffsetBuffer.data()+HaloOffsetBuffer.size()-Ngroups, sizeof(int), Ngroups, fd);
  // 	fseek(fd,sizeof(int)*Ngroups, SEEK_CUR);//skip offset
	  if(feof(fd))
	  {
		fprintf(stderr,"error:End-of-File in %s\n",filename);
		fflush(stderr);exit(1);  
	  }
	  fclose(fd);
	  return;
	}
 
  if(FileCounts>1)
	  sprintf(filename, filename_format.c_str(), "ids", iFile);
	else
	  sprintf(filename, filename_format.c_str(), "ids");
	
	myfopen(fd,filename,"r");
	int fileoffset;
	if(IsGroupV3)
	  fseek(fd, sizeof(Ngroups)+sizeof(TotNgroups)+sizeof(Nids)+sizeof(TotNids)+sizeof(NFiles)+sizeof(fileoffset), SEEK_CUR);
	else
	  fseek(fd, sizeof(Ngroups)+sizeof(Nids)+sizeof(TotNgroups)+sizeof(NFiles), SEEK_CUR);
	
	assert(end_particle<=Nids);
	if(end_particle<0) end_particle=Nids;
	HBTInt Nread=end_particle-start_particle;
	vector <PIDtype_t> PIDs(Nread);
	fseek(fd, sizeof(PIDtype_t)*start_particle, SEEK_CUR);
	myfread(PIDs.data(), sizeof(PIDtype_t), Nread, fd);
	fseek(fd, sizeof(PIDtype_t)*(Nids-end_particle), SEEK_CUR);
	if(long int extra_bytes=BytesToEOF(fd))
	{
	  cerr<<"Error: unexpected format of "<<filename<<endl;
	  cerr<<extra_bytes<<" extra bytes beyond particleId block! Check if particleId has a different type?\n";
	  exit(1);
	}
	if(feof(fd))
	{
	  fprintf(stderr,"error:End-of-File in %s\n",filename);
	  fflush(stderr);exit(1);  
	}
	fclose(fd);
  
  if(HBTConfig.ParticleIdRankStyle)
  {
	cerr<<"ParticleIdRankStyle not implemented for group files\n";
	exit(1);
  } else
  {
	HaloParticleBuffer.insert(HaloParticleBuffer.end(), PIDs.begin(), PIDs.end());
  }
}


#ifdef TEST_halo_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  mpi::environment env;
  mpi::communicator world;
  
  int snapshot_start, snapshot_end;
  if(0==world.rank())
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  broadcast(world, HBTConfig, 0);
  broadcast(world, snapshot_start, 0);
  broadcast(world, snapshot_end, 0);
  
  HaloSnapshot_t halo;
  halo.Load(world, HBTConfig.MaxSnapshotIndex);
  cout<<" Halo 0 from thread "<<world.rank()<<":"<<halo.Halos[0].Particles.size()<<", "<<halo.Halos[0].Particles[0]<<endl;
  return 0;
}
#endif
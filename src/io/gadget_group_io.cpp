#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <glob.h>
#include <climits>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstdio>
// #include <cmath>

#include "../mymath.h"
#include "../halo.h"
#include "gadget_group_io.h"

#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
#define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)

#define READ_META_DATA 0
#define READ_LEN_OFFSET 1
#define READ_PARTICLES 2
/*
inline bool IsNullParticle(const Particle_t &p)
{
  return p.Id==SpecialConst::NullParticleId;
}
  */
namespace GadgetGroup
{
bool GetGroupFileByteOrder(const char *filename, const int FileCounts)
{
  /* to check whether byteswap is needed, return true if yes, false if no, exit if error*/
  int Nfiles,n,ns;
  long offset;
  FILE *fp;
  string GroupFileFormat=HBTConfig.GroupFileFormat;
  
  if(GroupFileFormat=="gadget4")
	n=sizeof(GroupV4Header_t);
  else
	n=FileCounts;
  
  ns=n;
  swap_Nbyte(&ns,1,sizeof(ns));	
  
  myfopen(fp,filename,"r");
  
  if(GroupFileFormat=="gadget4")
	offset=0;
  else if(GroupFileFormat=="gadget3_int"||GroupFileFormat=="gadget3_long")
	offset=3*sizeof(int)+sizeof(long long);  
  else if(GroupFileFormat=="gadget3_MXXL")
    offset=2*sizeof(int)+2*sizeof(long long);
  else// if(GroupFileFormat=="gadget2_int"||GroupFileFormat=="gadget2_long")
	offset=3*sizeof(int); 

  fseek(fp,offset,SEEK_SET);
  size_t tmp_size=fread(&Nfiles,sizeof(int),1,fp);
  fclose(fp);
  
  if(Nfiles==n) return false;
  if(Nfiles==ns) return true;
  
  cerr<<"endianness check failed for: "<<filename<<", file format not expected:"<<Nfiles<<';'<<n<<';'<<ns<<endl;
  exit(1);
}

void GetFileNameFormat(int SnapshotId, string &format, int &FileCounts, bool &IsSubFile, bool &NeedByteSwap)
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
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts);
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
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts);
	format=fmt;
	return;
  }
  
  sprintf(basefmt, "%s/snapdir_%03d/group_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId,SnapshotId);
  sprintf(fmt, "%s.%%d", basefmt);
  sprintf(filename, fmt, "tab", ifile);
  if(file_exist(filename))
  {
	sprintf(pattern, basefmt, "tab");
	strcat(pattern, ".*");
	FileCounts=count_pattern_files(pattern);
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts);
	format=fmt;
	return;
  }
  
  sprintf(fmt, "%s/subhalo_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId);
  sprintf(filename, fmt, "tab");
  if(file_exist(filename))
  {
	IsSubFile=true;
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts);
	format=fmt;
	return;
  }
  
  sprintf(fmt, "%s/group_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId);
  sprintf(filename, fmt, "tab");
  if(file_exist(filename))
  {
	NeedByteSwap=GetGroupFileByteOrder(filename, FileCounts);
	format=fmt;
	return;
  }
  
  // 	return 0;
  cerr<<"Error: no group files found under "<<HBTConfig.HaloPath<<" at snapshot "<<SnapshotId<<endl;
  exit(1);
}  
  
class GroupFileReader_t
{
private:
  string filename_format;
  bool IsSubFile, NeedByteSwap;
  int iFile;
  int SnapshotId;
 
  template <class PIDtype_t>
  void ReadV2V3(int read_level, HBTInt start_particle, HBTInt end_particle);
  void ReadMXXL(int read_level, HBTInt start_particle, HBTInt end_particle);
  
public: 
  HBTInt NumberOfGroups;//in file, not necessarily the same as read
  HBTInt NumberOfParticles;//in file, not necessarily the same as read
  vector <HBTInt> Particles;
  vector <int> Len, Offset;
  int FileCounts;
  
  GroupFileReader_t(int SnapshotId): SnapshotId(SnapshotId), iFile(0), NumberOfGroups(0), NumberOfParticles(0), Particles(), Len(), Offset()
  {
	GetFileNameFormat(SnapshotId, filename_format, FileCounts, IsSubFile, NeedByteSwap);
  }
  void Read(int ifile, int read_level, HBTInt start_particle=0, HBTInt end_particle=-1);
};
void GroupFileReader_t::ReadMXXL(int read_level, HBTInt start_particle, HBTInt end_particle)
/* For reading MXXL (or PMO6K ) groups with compressed particle ids.
 *read_level: 0: meta data only ; 1: read group len and offset; 2: read particle list*/
{
  FILE *fd;
  char filename[1024];
  int Ngroups, Nids, NFiles;
  long long TotNids, TotNgroups; //long long TotNgroups, differ from V3
     
  if(HBTConfig.ParticleIdRankStyle)
  {
	cerr<<"ParticleIdRankStyle not implemented for group files\n";
	exit(1);
  } 
  
	if(FileCounts>1)
	  sprintf(filename, filename_format.c_str(), "tab", iFile);
	else
	  sprintf(filename, filename_format.c_str(), "tab");
		
	myfopen(fd,filename,"r");
	myfread(&Ngroups, sizeof(Ngroups), 1, fd);
    myfread(&TotNgroups, sizeof(TotNgroups), 1, fd);
    myfread(&Nids, sizeof(Nids), 1, fd);
    myfread(&TotNids,sizeof(TotNids),1,fd);
	
	NumberOfGroups=Ngroups;
	NumberOfParticles=Nids;
	myfread(&NFiles, sizeof(NFiles), 1, fd);
	if(IsSubFile)
	{
	  int Nsub;
      long long TotNsub;
	  myfread(&Nsub,sizeof(Nsub),1,fd);
	  myfread(&TotNsub,sizeof(TotNsub),1,fd);
	}
	if(FileCounts!=NFiles)
	{
	  cout<<"File count mismatch for file "<<filename<<": expect "<<FileCounts<<", got "<<NFiles<<endl;
	  exit(1);
	}
	if(read_level==READ_META_DATA) 
	{
	  fclose(fd);
	  return;
	}
	
	if(read_level==READ_LEN_OFFSET)
	{
	  Len.resize(Ngroups);
	  Offset.resize(Ngroups);
	  myfread(Len.data(), sizeof(int), Ngroups, fd);
	  myfread(Offset.data(), sizeof(int), Ngroups, fd);
  // 	fseek(fd,sizeof(int)*Ngroups, SEEK_CUR);//skip offset
	  if(feof(fd))
	  {
		fprintf(stderr,"error:End-of-File in %s\n",filename);
		fflush(stderr);exit(1);  
	  }
	  fclose(fd);
	  return;
	}
	fclose(fd);
 
  if(FileCounts>1)
	  sprintf(filename, filename_format.c_str(), "ids", iFile);
	else
	  sprintf(filename, filename_format.c_str(), "ids");
	
	myfopen(fd,filename,"r");
	long long fileoffset;//differ from V3
    fseek(fd, sizeof(Ngroups)+sizeof(TotNgroups)+sizeof(Nids)+sizeof(TotNids)+sizeof(NFiles)+sizeof(fileoffset), SEEK_CUR);
	
	assert(end_particle<=Nids);
	if(end_particle<0) end_particle=Nids;
	HBTInt Nread=end_particle-start_particle;
    
    typedef unsigned long long MyIDType;
    typedef MyIDType BucketType_t;
    #define BUCKET_SIZE (8*sizeof(BucketType_t))
    #define OBJECT_SIZE  40 //each 40bit id is an object  
    
    /* book-keeping for buckets */
    HBTInt start_bucket= (((long long) start_particle)* OBJECT_SIZE)/BUCKET_SIZE;
    HBTInt end_bucket= (((long long) end_particle)* OBJECT_SIZE)/BUCKET_SIZE;
    if((((long long) end_particle) * OBJECT_SIZE) % BUCKET_SIZE) 
        end_bucket++;
    int Nbucket = end_bucket-start_bucket;    
    int bucket_bits_start = (((long long) start_particle)* OBJECT_SIZE)%BUCKET_SIZE;
    int bucket_bits_left = BUCKET_SIZE-bucket_bits_start;

    int Nbucket_file;
    myfread(&Nbucket_file, sizeof(int), 1, fd);
    assert(Nbucket <= Nbucket_file);
    
    {
    /* read in the buckets containing the compressed ids */  
    vector <BucketType_t> comp(Nbucket);
    fseek(fd, sizeof(BucketType_t)*start_bucket, SEEK_CUR);
    myfread(comp.data(), sizeof(BucketType_t), Nbucket, fd);

    /* now do decompression */

    #define SET_MASK(m, n) ( ((((unsigned long long)1) << m)-1) << n ) //set m bits to 1 starting at offset n
    #define TAKE_BITS(x, m, n) ((((MyIDType) x) & SET_MASK(m, n)) >> n) //take m bits starting at bit n from data x
    
    Particles.resize(Nread);
    
    for(HBTInt i_object = 0, i_bucket=0; i_object < Nread; i_object++)
    {
        MyIDType ids_rec = 0;
        int object_bits_left = OBJECT_SIZE;
        int object_bits_start = 0;
        while(object_bits_left != 0)
        {
        int active_bits = min(bucket_bits_left, object_bits_left);
        //copy from bucket to object:
        ids_rec |= ( TAKE_BITS(comp[i_bucket], active_bits, bucket_bits_start) << object_bits_start);
        object_bits_left -= active_bits;
        object_bits_start += active_bits;
        bucket_bits_left -= active_bits;
        bucket_bits_start += active_bits;
        if(bucket_bits_left == 0)
            {
            i_bucket++;
            bucket_bits_left = BUCKET_SIZE;
            bucket_bits_start = 0;
            }
        }

        Particles[i_object] = ids_rec;
    }
    }
    
	fseek(fd, sizeof(BucketType_t)*(Nbucket_file-end_bucket), SEEK_CUR);
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
}

template <class PIDtype_t>
void GroupFileReader_t::ReadV2V3(int read_level, HBTInt start_particle, HBTInt end_particle)
/*read_level: 0: meta data only ; 1: read group len and offset; 2: read particle list*/
{
  FILE *fd;
  char filename[1024];
  int Ngroups, TotNgroups, Nids, NFiles;
  long long TotNids;
  bool IsGroupV3=("gadget3_int"==HBTConfig.GroupFileFormat||"gadget3_long"==HBTConfig.GroupFileFormat);
     
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
	NumberOfGroups=Ngroups;
	NumberOfParticles=Nids;
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
	if(read_level==READ_META_DATA) 
	{
	  fclose(fd);
	  return;
	}
	
	if(read_level==READ_LEN_OFFSET)
	{
	  Len.resize(Ngroups);
	  Offset.resize(Ngroups);
	  myfread(Len.data(), sizeof(int), Ngroups, fd);
	  myfread(Offset.data(), sizeof(int), Ngroups, fd);
  // 	fseek(fd,sizeof(int)*Ngroups, SEEK_CUR);//skip offset
	  if(feof(fd))
	  {
		fprintf(stderr,"error:End-of-File in %s\n",filename);
		fflush(stderr);exit(1);  
	  }
	  fclose(fd);
	  return;
	}
	fclose(fd);
 
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
	PIDtype_t Mask=HBTConfig.GroupParticleIdMask;
	if(Mask) 
	#pragma omp parallel for
	for(HBTInt i=0;i<Nread;i++)
	  PIDs[i]&=Mask;
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
	Particles.assign(PIDs.begin(), PIDs.end());
  }
}
void GroupFileReader_t::Read(int ifile, int read_level, HBTInt start_particle, HBTInt end_particle)
{
  iFile=ifile;
  
  string GroupFileFormat=HBTConfig.GroupFileFormat;
  if(GroupFileFormat=="gadget2_int"||GroupFileFormat=="gadget3_int")
      ReadV2V3<int>(read_level, start_particle, end_particle);
  else if(GroupFileFormat=="gadget2_long"||GroupFileFormat=="gadget3_long")
      ReadV2V3<long>(read_level, start_particle, end_particle);
  else if(GroupFileFormat=="gadget3_MXXL")
      ReadMXXL(read_level, start_particle, end_particle);
  else
      throw(runtime_error("unknown GroupFileFormat "+GroupFileFormat));
}

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
void Load(MpiWorker_t & world, int SnapshotId, vector <Halo_t> &Halos)
{
  GroupFileReader_t Reader(SnapshotId);
  int FileCounts=Reader.FileCounts;
  CountBuffer_t HaloLenBuffer;
  
  /*distribute tasks*/
  vector <FileAssignment_t> alltasks;
  FileAssignment_t thistask;
  if(world.rank()==0)
  {
	HBTInt Ngroups=0,Nparticles=0;
	CountBuffer_t HaloOffsetBuffer, FileOffset;
	FileOffset.resize(FileCounts);
	for(int iFile=0;iFile<FileCounts;iFile++)
	{
	  FileOffset[iFile]=Nparticles;
	  Reader.Read(iFile, READ_META_DATA);
	  Ngroups+=Reader.NumberOfGroups;
	  Nparticles+=Reader.NumberOfParticles;
	}
	HaloLenBuffer.reserve(Ngroups);
	HaloOffsetBuffer.reserve(Ngroups+1);
	for(int iFile=0;iFile<FileCounts;iFile++)
	{
	  Reader.Read(iFile, READ_LEN_OFFSET);
	  HaloLenBuffer.insert(HaloLenBuffer.end(), Reader.Len.begin(), Reader.Len.end());
	  HaloOffsetBuffer.insert(HaloOffsetBuffer.end(), Reader.Offset.begin(), Reader.Offset.end());
	}
	assert(CompileOffsets(HaloLenBuffer, HaloOffsetBuffer)==Nparticles);//the offsets in the group file may not always be correct. recompile.
	HaloOffsetBuffer.push_back(Nparticles);//end offset
	alltasks.resize(world.size());
	int nworkers=world.size();
	HBTInt npart_begin=0; 
	for(int rank=0;rank<world.size();rank++)
	  AssignHaloTasks(nworkers--, Nparticles, HaloOffsetBuffer, FileOffset, npart_begin, alltasks[rank]);
  }
  MPI_Datatype MPI_FileAssignment_t;
  FileAssignment_t().create_MPI_type(MPI_FileAssignment_t);
  MPI_Scatter(alltasks.data(), 1, MPI_FileAssignment_t, &thistask, 1, MPI_FileAssignment_t, 0, world.Communicator);
  MPI_Type_free(&MPI_FileAssignment_t);
  if(world.rank()==0)
  {
	for(int rank=1;rank<world.size();rank++)
	  MPI_Send(HaloLenBuffer.data()+alltasks[rank].haloid_begin, alltasks[rank].nhalo, MPI_HBT_INT, rank, 0, world.Communicator);
	HaloLenBuffer.resize(thistask.nhalo);
  }
  else
  {
	HaloLenBuffer.resize(thistask.nhalo);
	MPI_Recv(HaloLenBuffer.data(), HaloLenBuffer.size(), MPI_HBT_INT, 0, 0, world.Communicator, MPI_STATUS_IGNORE);
  }
  
  /* read particles*/
  ParticleIdBuffer_t ParticleBuffer;
  ParticleBuffer.reserve(thistask.npart);
  for(int i=0, ireader=0;i<world.size();i++, ireader++)
  {
	if(ireader==HBTConfig.MaxConcurrentIO) 
	{
	  ireader=0;//reset reader count
	  MPI_Barrier(world.Communicator);//wait for every thread to arrive.
	}
	if(i==world.rank())//read
	{
	  for(int iFile=thistask.ifile_begin;iFile<thistask.ifile_end;iFile++)
	  {
		HBTInt start_particle=0, end_particle=-1;
		if(iFile==thistask.ifile_begin) start_particle=thistask.firstfile_begin_part;
		if(iFile==thistask.ifile_end-1) end_particle=thistask.lastfile_end_part;
		Reader.Read(iFile, READ_PARTICLES, start_particle, end_particle);
		ParticleBuffer.insert(ParticleBuffer.end(), Reader.Particles.begin(), Reader.Particles.end());
	  }
	}
  }
  
  /* populate haloes*/
  Halos.resize(thistask.nhalo);
  auto p=ParticleBuffer.data();
  for(HBTInt i=0;i<Halos.size();i++)
  {
	Halos[i].HaloId=thistask.haloid_begin+i;
	Halos[i].Particles.resize(HaloLenBuffer[i]);
	for(HBTInt j=0;j<HaloLenBuffer[i];j++)
	  Halos[i].Particles[j].Id=*(p++);
  }
  
//   HBTInt np_node=0;
//   for(auto &&np: HaloLenBuffer)
// 	np_node+=np;//local
//   assert(np_node==thistask.npart);
  
  /*clear up buffers*/
  ParticleIdBuffer_t().swap(ParticleBuffer);
  CountBuffer_t().swap(HaloLenBuffer);
    
//   cout<<"Finished reading "<<Halos.size()<<" groups ("<<thistask.npart<<" particles) from file "<<thistask.ifile_begin<<" to "<<thistask.ifile_end-1<<" (total "<<FileCounts<<" files) on thread "<<world.rank()<<endl;  
//   if(Halos.size())
// 	cout<<"Halos loaded: "<<Halos.front().HaloId<<"-"<<Halos.back().HaloId<<endl; 
}

bool IsGadgetGroup(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 6)=="gadget";
}

}

#ifdef TEST_halo_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MpiWorker_t world(MPI_COMM_WORLD);
   
  int snapshot_start, snapshot_end;
  if(0==world.rank())
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  HBTConfig.BroadCast(world, 0, snapshot_start, snapshot_end);

  HaloSnapshot_t halo;
//   ParticleSnapshot_t snap;
//   snap.Load(world, snapshot_start);
//   cout<<"snapshot loaded\n";
  halo.Load(world, snapshot_start);
//   halo.UpdateParticles(world, snap);
  if(halo.Halos.size()>1)
  {
	auto & h=halo.Halos[1];
	cout<<" Halo 1 from thread "<<world.rank()<<":"<<"id="<<h.HaloId<<","<<h.Particles.size()<<", "<<h.Particles[5]<<endl;
  }
  
  MPI_Finalize();
  return 0;
}
#endif

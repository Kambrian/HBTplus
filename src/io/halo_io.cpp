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
  
  sprintf(basefmt, "%s/snapdir_%03d/group_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId,SnapshotId);
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

void HaloSnapshot_t::Load(int snapshot_index)
{
  SetSnapshotIndex(snapshot_index);
  switch(HBTConfig.GroupFileVariant)
  {
	case GROUP_FORMAT_GADGET2_INT:
	case GROUP_FORMAT_GADGET3_INT:
	  LoadGroupV2V3<int>();
	  break;
	case GROUP_FORMAT_GADGET2_LONG:
	case GROUP_FORMAT_GADGET3_LONG:
	  LoadGroupV2V3<long>();
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
void HaloSnapshot_t::LoadGroupV2V3()
{
  FILE *fd;
  char filename[1024];
  vector <int> Len;
  vector <HBTInt> Offset;
  int Ngroups, TotNgroups, Nids, NFiles;
  HBTInt NumberOfParticles,NumberOfHaloes;
  long long TotNids;
  bool NeedByteSwap;
  bool IsGroupV3=IS_GROUP_V3;
  
  string filename_format;
  bool IsSubFile;
  int FileCounts;
  GetFileNameFormat(filename_format, FileCounts, IsSubFile, NeedByteSwap);
  
  long long Nload=0;
NumberOfParticles=0;  
  for(int iFile=0;iFile<FileCounts;iFile++)
  {
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
	myfread(&NFiles, sizeof(NFiles), 1, fd);
	NumberOfParticles+=Nids;
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
	
	if(0==iFile)
	{
	  NumberOfHaloes=TotNgroups; 	  
	  Len.resize(TotNgroups);
	  Offset.resize(TotNgroups);
	}
	
	myfread(Len.data()+Nload, sizeof(int), Ngroups, fd);
	fseek(fd, sizeof(int)*Ngroups, SEEK_CUR);//skip offset
// 	myfread(Offset.data()+Nload, sizeof(int), Ngroups, fd);
	if(feof(fd))
	{
	  fprintf(stderr,"error:End-of-File in %s\n",filename);
	  fflush(stderr);exit(1);  
	}
	Nload+=Ngroups;
	fclose(fd);
  }
  if(Nload!=NumberOfHaloes)
  {
	cerr<<"error:Num groups loaded not match: "<<Nload<<','<<NumberOfHaloes<<"\n";
	exit(1);
  }
  assert(CompileOffsets(Len, Offset)==NumberOfParticles);//compile offsets in case the stored offsets are wrong.
  cout<<"GroupSnap="<<SnapshotIndex<<" ("<<SnapshotId<<") TotNids="<<NumberOfParticles<<" TotNgroups="<<TotNgroups<<" NFiles="<<NFiles<<endl;
  
  vector <PIDtype_t> PIDs(NumberOfParticles);
  
  Nload=0;  
  for(int iFile=0;iFile<FileCounts;iFile++)
  {
	if(FileCounts>1)
	  sprintf(filename, filename_format.c_str(), "ids", iFile);
	else
	  sprintf(filename, filename_format.c_str(), "ids");
	
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
	myfread(&NFiles, sizeof(NFiles), 1, fd);
	//file offset:
	if(IsGroupV3)
	{
	  int file_offset;
	  myfread(&file_offset,sizeof(int),1,fd);
	}
	
	myfread(&PIDs[Nload], sizeof(PIDtype_t), Nids, fd);
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
	Nload+=Nids;
	fclose(fd);
  }
 		
  if(Nload!=NumberOfParticles)
  {
	cerr<<"error:Number of  group particles loaded not match: "<<Nload<<','<<NumberOfParticles<<endl;
	exit(1);
  }
  
  PIDtype_t Mask=HBTConfig.GroupParticleIdMask;
  if(Mask)
  #pragma omp parallel for
  for(HBTInt i=0;i<NumberOfParticles;i++)
	PIDs[i]&=Mask;
  
  TotNumberOfParticles=NumberOfParticles;
  Halos.resize(NumberOfHaloes);
  if(HBTConfig.ParticleIdRankStyle)
  {
	cerr<<"ParticleIdRankStyle not implemented for group files\n";
	exit(1);
  } else
  {
#pragma omp parallel for
	for(int i=0;i<NumberOfHaloes;i++)
		Halos[i].Particles.assign(PIDs.begin()+Offset[i], PIDs.begin()+Offset[i]+Len[i]);
// 	if(typeid(HBTInt)==typeid(PIDtype_t))//consider optimizing by using memcpy
  }
  
  NumPartOfLargestHalo=0;
  if(Len.size()) NumPartOfLargestHalo=*max_element(Len.begin(), Len.end());
  
  for(int i=1;i<NumberOfHaloes;i++)
  {
	if(Len[i]>Len[i-1])
	{
	  fprintf(stderr, "warning: groups not sorted with mass\n");
	  break;
	}
  }
}


#ifdef TEST_halo_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  HBTConfig.ParseConfigFile(argv[1]);
  HaloSnapshot_t halo;
  halo.Load(HBTConfig.MaxSnapshotIndex);
  cout<<halo.Halos.size()<<";"<<halo.TotNumberOfParticles<<endl;
  cout<<halo.Halos[10].Particles.size()<<endl;
  cout<<halo.Halos[10].Particles[0]<<endl;
/*  ParticleSnapshot_t snap;
  snap.Load(HBTConfig.MaxSnapshotIndex);
  halo.ParticleIdToIndex(snap);
  int ids[]={1,53};
  for(auto &hid: ids)
  {
	auto & h=halo.Halos[hid];
	auto pid=h.Particles[5];
	const HBTxyz &pos=snap.GetComovingPosition(pid);
	const HBTxyz &vel=snap.GetPhysicalVelocity(pid);
	cout<<" Halo id="<<hid<<","<<h.Particles.size()<<", ["<<snap.GetParticleId(pid)<<","<<snap.GetParticleMass(pid)<<",("<<pos[0]<<","<<pos[1]<<","<<pos[2]<<"), ("<<vel[0]<<","<<vel[1]<<","<<vel[2]<<") ]"<<endl;
  }
  */
  return 0;
}
#endif
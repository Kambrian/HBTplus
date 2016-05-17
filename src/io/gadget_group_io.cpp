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
#include "gadget_group_io.h"

#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
#define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)

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

template <class PIDtype_t>
HBTInt LoadGroupV2V3(int SnapshotId, vector <Halo_t> &Halos)
{
  FILE *fd;
  char filename[1024];
  vector <int> Len;
  vector <HBTInt> Offset;
  int Ngroups, TotNgroups, Nids, NFiles;
  HBTInt NumberOfParticles,NumberOfHaloes;
  long long TotNids;
  bool NeedByteSwap;
  bool IsGroupV3=("gadget3_int"==HBTConfig.GroupFileFormat||"gadget3_long"==HBTConfig.GroupFileFormat);
  
  string filename_format;
  bool IsSubFile;
  int FileCounts;
  GetFileNameFormat(SnapshotId, filename_format, FileCounts, IsSubFile, NeedByteSwap);
  
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
  cout<<"GroupSnap="<<SnapshotId<<" TotNids="<<NumberOfParticles<<" TotNgroups="<<TotNgroups<<" NFiles="<<NFiles<<endl;
  
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
  
  for(int i=1;i<NumberOfHaloes;i++)
  {
	if(Len[i]>Len[i-1])
	{
	  fprintf(stderr, "warning: groups not sorted with mass\n");
	  break;
	}
  }
  
  return NumberOfParticles;
}

HBTInt Load(int SnapshotId, vector <Halo_t> &Halos)
{
  string GroupFileFormat=HBTConfig.GroupFileFormat;
  HBTInt np=0;
  if(GroupFileFormat=="gadget2_int"||GroupFileFormat=="gadget3_int")
	np=GadgetGroup::LoadGroupV2V3<int>(SnapshotId, Halos);
  else if(GroupFileFormat=="gadget2_long"||GroupFileFormat=="gadget3_long")
	np=GadgetGroup::LoadGroupV2V3<long>(SnapshotId, Halos);
  else
	throw(runtime_error("unknown GroupFileFormat "+GroupFileFormat));
  
  return np;
}

bool IsGadgetGroup(const string &GroupFileFormat)
{
  return GroupFileFormat.substr(0, 6)=="gadget";
}

}
using namespace std;
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

#include "../mymath.h"
#include "halo.h"

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
	  offset=5*4;  //3*int+1*longlong
	  break;
	default:
	  offset=3*4;  //3*int
  }  
  fseek(fp,offset,SEEK_SET);
  fread(&Nfiles,sizeof(int),1,fp);
  fclose(fp);
  
  if(Nfiles==n) return false;
  if(Nfiles==ns) return true;
  
  cerr<<"endianness check failed for: "<<filename<<", file format not expected:"<<Nfiles<<';'<<n<<';'<<ns<<endl;
  exit(1);
}

void HaloSnapshot_t::GetFileNameFormat(Parameter_t &param, string &format, int &FileCounts, bool &IsSubFile, bool &NeedByteSwap)
{
  IsSubFile=false;
  FileCounts=1;
  char buf[1024], pattern[1024], basefmt[1024], fmt[1024];
  const int ifile=0;
  sprintf(basefmt, "%s/groups_%03d/subhalo_%%s_%03d",param.HaloPath.c_str(),SnapshotId,SnapshotId);
  sprintf(fmt, "%s.%%d", basefmt);
  sprintf(buf, fmt, "tab", ifile);
  if(file_exist(buf))
  {
	IsSubFile=true;
	sprintf(pattern, "%s.*", basefmt);
	sprintf(pattern, pattern, "tab");
	FileCounts=count_pattern_files(pattern);
	NeedByteSwap=GetGroupFileByteOrder(buf, FileCounts, param.GroupFileVariant);
	return;
  }
  
  sprintf(basefmt, "%s/groups_%03d/group_%%s_%03d",param.HaloPath.c_str(),SnapshotId,SnapshotId);
  sprintf(fmt, "%s.%%d", basefmt);
  sprintf(buf, fmt, "tab", ifile);
  if(file_exist(buf))
  {
	sprintf(pattern, "%s.*", basefmt);
	sprintf(pattern, pattern, "tab");
	FileCounts=count_pattern_files(pattern);
	NeedByteSwap=GetGroupFileByteOrder(buf, FileCounts, param.GroupFileVariant);
	return;
  }
  
  sprintf(fmt, "%s/subhalo_%%s_%03d",param.HaloPath.c_str(),SnapshotId);
  sprintf(buf, fmt, "tab");
  if(file_exist(buf))
  {
	IsSubFile=true;
	NeedByteSwap=GetGroupFileByteOrder(buf, FileCounts, param.GroupFileVariant);
	return;
  }
  
  sprintf(fmt, "%s/group_%%s_%03d",param.HaloPath.c_str(),SnapshotId);
  sprintf(buf, fmt, "tab");
  if(file_exist(buf))
  {
	NeedByteSwap=GetGroupFileByteOrder(buf, FileCounts, param.GroupFileVariant);
	return;
  }
  
  // 	return 0;
  cerr<<"Error: no group files found under "<<param.HaloPath<<" at snapshot "<<SnapshotId<<endl;
  exit(1);
}

void HaloSnapshot_t::Load(Parameter_t &param, int snapshot_index)
{
  SetSnapshotIndex(param, snapshot_index);
  switch(param.GroupFileVariant)
  {
	case GROUP_FORMAT_GADGET3_INT:
	  int i;
	  LoadGroupV3(param, i);
	  break;
	case GROUP_FORMAT_GADGET3_LONG:
	  long j;
	  LoadGroupV3(param, j);
	  break;
	default:
	  cerr<<"GroupFileVariant="<<param.GroupFileVariant<<" not implemented yet.\n";
	  exit(1);
  }
}

void HaloSnapshot_t::Clear()
{
  delete [] ParticleList;
  delete [] HaloLength;
  delete [] HaloOffset;
  NumberOfHaloes=0;
  NumberOfParticles=0;
}

template <class PIDtype_t>
void HaloSnapshot_t::LoadGroupV3(Parameter_t &param, PIDtype_t dummy)
{
  FILE *fd;
  char buf[1024];
  int *Len, *Offset;
  int Ngroups, TotNgroups, Nids, NFiles;
  long long TotNids;
  bool NeedByteSwap;
  
  string filename_format;
  bool IsSubFile;
  int FileCounts;
  GetFileNameFormat(param, filename_format, FileCounts, IsSubFile, NeedByteSwap);
  
  long long Nload=0;  
  for(int iFile=0;iFile<FileCounts;iFile++)
  {
	if(FileCounts>1)
	  sprintf(buf, filename_format.c_str(), "tab", iFile);
	else
	  sprintf(buf, filename_format.c_str(), "tab");
		
	myfopen(fd,buf,"r");
	myfread(&Ngroups, sizeof(Ngroups), 1, fd);
	myfread(&TotNgroups, sizeof(TotNgroups), 1, fd);
	myfread(&Nids, sizeof(Nids), 1, fd);
	myfread(&TotNids,sizeof(TotNids),1,fd);
	myfread(&NFiles, sizeof(NFiles), 1, fd);
	if(IsSubFile)
	{
	  int Nsub,TotNsub;
	  myfread(&Nsub,sizeof(int),1,fd);
	  myfread(&TotNsub,sizeof(int),1,fd);
	}
	if(FileCounts!=NFiles)
	{
	  fprintf(stderr,"error: number of grpfiles specified not the same as stored: %d,%d\n for file %s\n",
			  FileCounts,NFiles,buf);
	  fflush(stderr);
	  exit(1);
	}
	
	if(0==iFile)
	{
	  NumberOfHaloes=TotNgroups;
	  NumberOfParticles=TotNids;
	  if(sizeof(HBTInt)==4 && TotNids>INT_MAX)
	  {
		fprintf(stderr,"error: TotNids larger than HBTInt can hold! %lld\n",TotNids);
		exit(1);
	  }
	  fprintf(stdout,"Snap=%d (%d)  TotNids=%ld  TotNgroups=%d  NFiles=%d\n", SnapshotIndex, SnapshotId, (long)NumberOfParticles, (int)NumberOfHaloes, NFiles);
	  
	  Len=new int[NumberOfHaloes];
	  Offset=new int[NumberOfHaloes];
	}
	
	myfread(Len+Nload, sizeof(int), Ngroups, fd);
	myfread(Offset+Nload, sizeof(int), Ngroups, fd);
	if(feof(fd))
	{
	  fprintf(stderr,"error:End-of-File in %s\n",buf);
	  fflush(stderr);exit(1);  
	}
	Nload+=Ngroups;
	fclose(fd);
  }
  if(Nload!=NumberOfHaloes)
  {
	fprintf(stderr,"error:Num groups loaded not match: %lld,"HBTIFMT"\n",Nload,NumberOfHaloes);
	fflush(stderr);exit(1);
  }
  if(typeid(HBTInt)==typeid(int))
  {
	HaloLength=(HBTInt *)Len;
	HaloOffset=(HBTInt *)Offset;
  }
  else
  {
	HaloLength=new HBTInt[NumberOfHaloes];
	HaloOffset=new HBTInt[NumberOfHaloes];
	for(int i=0;i<NumberOfHaloes;i++)
	{
	  HaloLength[i]=Len[i];
	  HaloOffset[i]=Offset[i];
	}
	delete [] Len;
	delete [] Offset;
  }
  
  PIDtype_t *PID=new PIDtype_t[NumberOfParticles];
  
  Nload=0;  
  for(int iFile=0;iFile<FileCounts;iFile++)
  {
	if(FileCounts>1)
	  sprintf(buf, filename_format.c_str(), "ids", iFile);
	else
	  sprintf(buf, filename_format.c_str(), "ids");
	
	myfopen(fd,buf,"r");
	myfread(&Ngroups, sizeof(Ngroups), 1, fd);
	myfread(&TotNgroups, sizeof(TotNgroups), 1, fd);
	myfread(&Nids, sizeof(Nids), 1, fd);
	myfread(&TotNids,sizeof(TotNids),1,fd);
	myfread(&NFiles, sizeof(NFiles), 1, fd);
	//file offset:
	int dummy;
	myfread(&dummy,sizeof(int),1,fd);
	
	myfread(PID+Nload, sizeof(PIDtype_t), Nids, fd);
	if(feof(fd))
	{
	  fprintf(stderr,"error:End-of-File in %s\n",buf);
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
  
  if(param.ParticleIdRankStyle)
  {
	cerr<<"ParticleIdRankStyle not implemented for group files\n";
	exit(1);
  }else
  {
	if(typeid(HBTInt)==typeid(PIDtype_t))
	  ParticleList=(HBTInt *)PID;
	else
	{
	  ParticleList=new HBTInt[NumberOfParticles];
	  for(HBTInt i=0;i<NumberOfParticles;i++)
		  ParticleList[i]=PID[i];
	  delete [] PID;
	}
  }
  for(int i=1;i<NumberOfHaloes;i++)
  {
	if(HaloLength[i]>HaloLength[i-1])
	{
	  fprintf(stderr, "warning: groups not sorted with mass\n");
	  break;
	}
  }
}	
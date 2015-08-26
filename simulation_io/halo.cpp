using namespace std;
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "../mymath.h"
#include "halo.h"

#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
#define ReadBlockSize(a) myfread(&a,sizeof(a),1,fp)

inline size_t SkipBlock(FILE *fp, bool NeedByteSwap)
{
  int blocksize,blocksize2;
  
  ReadBlockSize(blocksize);	  
  fseek(fp, blocksize, SEEK_CUR);
  ReadBlockSize(blocksize2);
  assert(blocksize==blocksize2);
  return blocksize;
}

void Halo_t::GetFileName(Parameter_t &param, int ifile, string &filename)
{
  FILE *fp;
  char buf[1024];
  
  sprintf(buf,"%s/snapdir_%03d/%s_%03d.%d",param.SnapshotPath.c_str(),SnapshotId,param.SnapshotFileBase.c_str(),SnapshotId,ifile);
  if(ifile==0)
	if(!file_exist(buf)) sprintf(buf,"%s/%s_%03d",param.SnapshotPath.c_str(),param.SnapshotFileBase.c_str(),SnapshotId); //try the other convention
  if(!file_exist(buf)) sprintf(buf,"%s/%s_%03d.%d",param.SnapshotPath.c_str(),param.SnapshotFileBase.c_str(),SnapshotId,ifile); //try the other convention
  if(!file_exist(buf)) sprintf(buf,"%s/%d/%s.%d",param.SnapshotPath.c_str(),SnapshotId,param.SnapshotFileBase.c_str(),ifile);//for BJL's RAMSES output
  if(!file_exist(buf))
  {
	cerr<<"Failed to find a snapshot file at "<<buf<<endl;
	exit(1);
  }
  filename=buf;
}
int find_group_file(HBTInt Nsnap, char *GrpPath, int *grpfile_type)
{
  int i=0;
  char buf[1024],pattern[1024];
    
  sprintf(buf, "%s/groups_%03d/subhalo_ids_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
  if(try_readfile(buf))
  {
	*grpfile_type=1;
	sprintf(pattern, "%s/groups_%03d/subhalo_tab_%03d.*",GrpPath,(int)Nsnap,(int)Nsnap);
	return count_pattern_files(pattern);
  }
  
  sprintf(buf, "%s/groups_%03d/group_tab_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
  if(try_readfile(buf))
  {
	*grpfile_type=2;
	sprintf(pattern, "%s/groups_%03d/group_tab_%03d.*",GrpPath,(int)Nsnap,(int)Nsnap);
	return count_pattern_files(pattern);
  }
  
  sprintf(buf, "%s/subhalo_tab_%03d",GrpPath,(int)Nsnap);
  if(try_readfile(buf))
  {
	  *grpfile_type=3;
	  return 1;
  }
  
  sprintf(buf, "%s/group_tab_%03d",GrpPath,(int)Nsnap);
  if(try_readfile(buf))
  {
	*grpfile_type=4;
	return 1;
  }
  
  return 0; //0 files matching
}
void load_group_catalogue_v3(HBTInt Nsnap,CATALOGUE *Cat,char *GrpPath)
{//PGADGET-3's subfind format
  FILE *fd;
  char buf[1024];
  int Ngroups,TotNgroups,Nids,NFiles, FileCounts;
  long long i,Nload,TotNids;
  int ByteOrder,grpfile_type;
  
struct 
{
int *Len;
int *Offset;
IDatInt *PID; 
}ICat;
  
  	#ifdef SNAPLIST
	Nsnap=snaplist[Nsnap];
	#endif

  Nload=0;
  FileCounts=find_group_file(Nsnap, GrpPath, &grpfile_type);
  if(!FileCounts) 
  {
	fprintf(logfile, "Error: no group files found under %s at snapshot %d\n", GrpPath, Nsnap);
	exit(1);
  }
  
  for(i=0;i<FileCounts;i++)
  {
	switch(grpfile_type)
	{
	  case 1:
		sprintf(buf, "%s/groups_%03d/subhalo_tab_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
		break;
	  case 2:
		sprintf(buf, "%s/groups_%03d/group_tab_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
		break;
	  case 3:
		sprintf(buf, "%s/subhalo_tab_%03d",GrpPath,(int)Nsnap);
		break;
	  case 4:
		sprintf(buf, "%s/group_tab_%03d",GrpPath,(int)Nsnap);
		break;
	  default:
		fprintf(logfile, "Error: wrong grpfile_type=%d at snap %d\n", grpfile_type, Nsnap);
		exit(1);
	}
	
  ByteOrder=check_grpcat_byteorder(buf, FileCounts);
  	
  myfopen(fd,buf,"r");

  myfread(&Ngroups, sizeof(int), 1, fd);
  myfread(&TotNgroups, sizeof(int), 1, fd);
  Cat->Ngroups=TotNgroups;
  myfread(&Nids, sizeof(int), 1, fd);
  myfread(&TotNids,sizeof(long long),1,fd);
  #ifndef HBT_INT8
  if(TotNids>INT_MAX)
  {
	  printf("error: TotNids larger than HBTInt can hold! %lld\n",TotNids);
	  exit(1);
  }
  #endif
  Cat->Nids=TotNids;
  myfread(&NFiles, sizeof(int), 1, fd);
 
  if(1==grpfile_type||3==grpfile_type)
  {
  int Nsub,TotNsub;
  myfread(&Nsub,sizeof(int),1,fd);
  myfread(&TotNsub,sizeof(int),1,fd);
  }
  
  if(FileCounts!=NFiles)
	  {
		  fprintf(logfile,"error: number of grpfiles specified not the same as stored: %d,%d\n for file %s\n",
		  FileCounts,(int)NFiles,buf);
		  fflush(logfile);
		  exit(1);
	  }
 
  if(0==i)
  fprintf(logfile,"Snap="HBTIFMT"  Ngroups=%d  Nids=%d  TotNgroups=%d  NFiles=%d\n", Nsnap, Ngroups, Nids, (int)(Cat->Ngroups), NFiles);
  
  if(0==i)
  {
  ICat.Len= mymalloc(sizeof(int)*Cat->Ngroups);
  ICat.Offset=mymalloc(sizeof(int)*Cat->Ngroups);
  }
  
  myfread(ICat.Len+Nload, sizeof(int), Ngroups, fd);
  myfread(ICat.Offset+Nload, sizeof(int), Ngroups, fd);
  if(feof(fd))
  {
	fprintf(logfile,"error:End-of-File in %s\n",buf);
	fflush(logfile);exit(1);  
  }
  Nload+=Ngroups;
  fclose(fd);
  }
  if(Nload!=Cat->Ngroups)
  {
	  fprintf(logfile,"error:Num groups loaded not match: %lld,"HBTIFMT"\n",Nload,Cat->Ngroups);
	  fflush(logfile);exit(1);
  }
  
  ICat.PID=mymalloc(sizeof(IDatInt)*Cat->Nids);
  Cat->ID2Halo=mymalloc(sizeof(HBTInt)*NP_DM);//consider move this out.............
  Nload=0;
  for(i=0;i<FileCounts;i++)
  {
  switch(grpfile_type)
  {
  case 1:
	sprintf(buf, "%s/groups_%03d/subhalo_ids_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
	break;
  case 2:
	sprintf(buf, "%s/groups_%03d/group_ids_%03d.%d",GrpPath,(int)Nsnap,(int)Nsnap,(int)i);
	break;
  case 3:
	sprintf(buf, "%s/subhalo_ids_%03d",GrpPath,(int)Nsnap);
	break;
  case 4:
	sprintf(buf, "%s/group_ids_%03d",GrpPath,(int)Nsnap);
  break;
  default:
	fprintf(logfile,"error: grpfile_type not assigned? %s\n",buf);
	exit(1);
  }
  ByteOrder=check_grpcat_byteorder(buf, FileCounts);
  myfopen(fd,buf,"r");

  myfread(&Ngroups, sizeof(int), 1, fd);
  myfread(&TotNgroups, sizeof(int), 1, fd);
  Cat->Ngroups=TotNgroups;
  myfread(&Nids, sizeof(int), 1, fd);
  myfread(&TotNids,sizeof(long long),1,fd);
  #ifndef HBT_INT8
  if(TotNids>INT_MAX)
  {
	  printf("error: TotNids larger than HBTInt can hold! %lld\n",TotNids);
	  exit(1);
  }
  #endif
  Cat->Nids=TotNids;
  myfread(&NFiles, sizeof(int), 1, fd);
  //file offset:
  int dummy;
  myfread(&dummy,sizeof(int),1,fd);

  
  myfread(ICat.PID+Nload, sizeof(IDatInt), Nids, fd);
  if(feof(fd))
  {
	fprintf(logfile,"error:End-of-File in %s\n",buf);
	fflush(logfile);exit(1);  
  }
  Nload+=Nids;
  fclose(fd);
	}
	
  if(Nload!=Cat->Nids)
  {
	  fprintf(logfile,"error:Num grpparticles loaded not match: %lld,%lld\n",Nload,(long long)Cat->Nids);
	  fflush(logfile);
	  exit(1);
  }
  
  
#ifndef HBT_INT8
	  Cat->Len=ICat.Len;
	  Cat->Offset=ICat.Offset;
#else
	  Cat->Len=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
	  Cat->Offset=mymalloc(sizeof(HBTInt)*Cat->Ngroups);
	  for(i=0;i<Cat->Ngroups;i++)
	  {
		  Cat->Len[i]=ICat.Len[i];
		  Cat->Offset[i]=ICat.Offset[i];
	  }
	  myfree(ICat.Len);
	  myfree(ICat.Offset);
#endif
  
  #ifdef HBTPID_RANKSTYLE
  IDatInt *PIDs,*p;  
  PIDs=load_PIDs_Sorted();
  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
  for(i=0;i<Cat->Nids;i++)
  {	  
	p=bsearch(&(ICat.PID[i]),PIDs,NP_DM,sizeof(IDatInt),comp_IDatInt);
	Cat->PIDorIndex[i]=p-PIDs;
  }
  myfree(PIDs);
  myfree(ICat.PID);
  #else
#ifdef SAME_INTTYPE
	  Cat->PIDorIndex=ICat.PID;
#else
	  Cat->PIDorIndex=mymalloc(sizeof(HBTInt)*Cat->Nids);
	  for(i=0;i<Cat->Nids;i++)
	  Cat->PIDorIndex[i]=ICat.PID[i];
	  myfree(ICat.PID);
#endif 
  #endif
  
	for(i=1;i<Cat->Ngroups;i++)
	{
		if(Cat->Len[i]>Cat->Len[i-1])
		{
		printf("warning: groups not sorted with mass\n");
		}
	}
	
}	
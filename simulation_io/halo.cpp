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
	  offset=3*sizeof(int)+sizeof(long long);  
	  break;
	default:
	  offset=3*sizeof(int); 
  }  
  fseek(fp,offset,SEEK_SET);
  fread(&Nfiles,sizeof(int),1,fp);
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
  char buf[1024], pattern[1024], basefmt[1024], fmt[1024];
  const int ifile=0;
  sprintf(basefmt, "%s/groups_%03d/subhalo_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId,SnapshotId);
  sprintf(fmt, "%s.%%d", basefmt);
  sprintf(buf, fmt, "tab", ifile);
  if(file_exist(buf))
  {
	IsSubFile=true;
	sprintf(pattern, basefmt, "tab");
	strcat(pattern, ".*");
	FileCounts=count_pattern_files(pattern);
	NeedByteSwap=GetGroupFileByteOrder(buf, FileCounts, HBTConfig.GroupFileVariant);
	format=fmt;
	return;
  }
  
  sprintf(basefmt, "%s/groups_%03d/group_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId,SnapshotId);
  sprintf(fmt, "%s.%%d", basefmt);
  sprintf(buf, fmt, "tab", ifile);
  if(file_exist(buf))
  {
	sprintf(pattern, basefmt, "tab");
	strcat(pattern, ".*");
	FileCounts=count_pattern_files(pattern);
	NeedByteSwap=GetGroupFileByteOrder(buf, FileCounts, HBTConfig.GroupFileVariant);
	format=fmt;
	return;
  }
  
  sprintf(fmt, "%s/subhalo_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId);
  sprintf(buf, fmt, "tab");
  if(file_exist(buf))
  {
	IsSubFile=true;
	NeedByteSwap=GetGroupFileByteOrder(buf, FileCounts, HBTConfig.GroupFileVariant);
	format=fmt;
	return;
  }
  
  sprintf(fmt, "%s/group_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId);
  sprintf(buf, fmt, "tab");
  if(file_exist(buf))
  {
	NeedByteSwap=GetGroupFileByteOrder(buf, FileCounts, HBTConfig.GroupFileVariant);
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
  int dummy_integer=0;
  long dummy_long=0;
  switch(HBTConfig.GroupFileVariant)
  {
	case GROUP_FORMAT_GADGET3_INT:
	  LoadGroupV3(dummy_integer);
	  break;
	case GROUP_FORMAT_GADGET3_LONG:
	  LoadGroupV3(dummy_long);
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

template <class PIDtype_t>
void HaloSnapshot_t::LoadGroupV3(PIDtype_t dummy)
//the value of dummy is irrelevant, only its type is used.
{
  FILE *fd;
  char buf[1024];
  int *Len, *Offset;
  int Ngroups, TotNgroups, Nids, NFiles;
  HBTInt NumberOfParticles,NumberOfHaloes;
  long long TotNids;
  bool NeedByteSwap;
  
  string filename_format;
  bool IsSubFile;
  int FileCounts;
  GetFileNameFormat(filename_format, FileCounts, IsSubFile, NeedByteSwap);
  
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
	  fprintf(stdout,"Snap=%d (%d)  TotNids=%lld  TotNgroups=%d  NFiles=%d\n", SnapshotIndex, SnapshotId, TotNids, TotNgroups, NFiles);
	  
	  Len=new int[TotNgroups];
	  Offset=new int[TotNgroups];
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
	cerr<<"error:Num groups loaded not match: "<<Nload<<','<<NumberOfHaloes<<"\n";
	exit(1);
  }
  
  vector <PIDtype_t> PIDs(NumberOfParticles);
  
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
	
	myfread(&PIDs[Nload], sizeof(PIDtype_t), Nids, fd);
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
  
  TotNumberOfParticles=NumberOfParticles;
  Halos.resize(NumberOfHaloes);
  if(HBTConfig.ParticleIdRankStyle)
  {
	cerr<<"ParticleIdRankStyle not implemented for group files\n";
	exit(1);
  } else
  {
	for(int i=0;i<NumberOfHaloes;i++)
		Halos[i].Particles.assign(PIDs.begin()+Offset[i], PIDs.begin()+Offset[i]+Len[i]);
// 	if(typeid(HBTInt)==typeid(PIDtype_t))//consider optimizing by using memcpy
  }
  
  NumPartOfLargestHalo=*max_element(Len, Len+NumberOfHaloes);
  
  for(int i=1;i<NumberOfHaloes;i++)
  {
	if(Len[i]>Len[i-1])
	{
	  fprintf(stderr, "warning: groups not sorted with mass\n");
	  break;
	}
  }
}

void HaloSnapshot_t::ParticleIdToIndex(const Snapshot_t& snapshot)
{
  chrono::high_resolution_clock::time_point time_begin, time_end;
#pragma omp master
  {
  SnapshotPointer=&snapshot;
  time_begin = chrono::high_resolution_clock::now();
  }
#pragma omp for //maybe try collapse(2)? need to remove intermediate variables to enable this.
  for(HBTInt haloid=0;haloid<Halos.size();haloid++)
  {
	Halo_t::ParticleList_t & Particles=Halos[haloid].Particles;
	HBTInt np=Particles.size();
	for(HBTInt pid=0;pid<np;pid++)
	  Particles[pid]=snapshot.GetParticleIndex(Particles[pid]);//this should be safe since its a const func
  }
  #pragma omp master
  {
	time_end = chrono::high_resolution_clock::now();
	auto elapsed = chrono::duration_cast<chrono::duration<double>>(time_end - time_begin);
	cout << "Halo ParticleId lookup took " << elapsed.count() <<"seconds"<< endl;
  }
}

void HaloSnapshot_t::ParticleIndexToId()
{
#pragma omp for
  for(HBTInt haloid=0;haloid<Halos.size();haloid++)
  {
	Halo_t::ParticleList_t &Particles=Halos[haloid].Particles;
	HBTInt nP=Halos[haloid].Particles.size();
	for(HBTInt pid=0;pid<nP;pid++)
	  Particles[pid]=SnapshotPointer->GetParticleId(Particles[pid]);
  }
#pragma omp single
  SnapshotPointer=nullptr;
}

void HaloSnapshot_t::AverageCoordinates()
{
  for(HBTInt i=0;i<Halos.size();i++)
  {
	SnapshotPointer->AveragePosition(Halos[i].ComovingPosition, Halos[i].Particles.data(), Halos[i].Particles.size());
	SnapshotPointer->AverageVelocity(Halos[i].PhysicalVelocity, Halos[i].Particles.data(), Halos[i].Particles.size());
  }
}

#ifdef TEST_SIMU_IO_halo
#include "../config_parser.h"

int main(int argc, char **argv)
{
  HBTConfig.ParseConfigFile(argv[1]);
  HaloSnapshot_t halo;
  halo.Load(HBTConfig.MaxSnapshotIndex);
  cout<<halo.Halos.size()<<";"<<halo.TotNumberOfParticles<<endl;
  cout<<halo.Halos[10].Particles.size()<<endl;
  cout<<halo.Halos[10].Particles[0]<<endl;
  return 0;
}
#endif
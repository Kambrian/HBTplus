/*utility function for reading gadget output of halo virial mass and radius. 
 *This file is not used by main program of HBT+, and is provided for convenience of analysing gadget group data*/
#ifndef GADGET_VIRIAL_IO_HEADER_INCLUDED
#define GADGET_VIRIAL_IO_HEADER_INCLUDED

#include "gadget_group_io.h"
#define myfread(buf,size,count,fp) fread_swap(buf,size,count,fp,NeedByteSwap)
namespace GadgetVirialIO{
template <class MyReal>
void LoadVirial(int SnapshotId, vector <MyReal> &Mvir, vector <MyReal> &Rvir, const string & virtype)
{
//   typedef float MyReal;
  
  int itype;
  if(virtype=="Mean200")
    itype=0;
  else if(virtype=="Crit200")
    itype=1;
  else if(virtype=="TopHat")
    itype=2;
  else
    throw(runtime_error("unknow virtype"+virtype));
	
  FILE *fd;
  char filename[1024];
  vector <int> Len;
  vector <HBTInt> Offset;
  int Ngroups, TotNgroups, Nids, NFiles;
  HBTInt NumberOfHaloes;
  long long TotNids;
  bool NeedByteSwap;
  bool IsGroupV3=("gadget3_int"==HBTConfig.GroupFileFormat||"gadget3_long"==HBTConfig.GroupFileFormat);
  
  string filename_format;
  bool IsSubFile;
  int FileCounts;
  GadgetGroup::GetFileNameFormat(SnapshotId, filename_format, FileCounts, IsSubFile, NeedByteSwap);
  assert(IsSubFile);
  
  long long Nload=0;
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
	int Nsub,TotNsub;
	myfread(&Nsub,sizeof(int),1,fd);
	myfread(&TotNsub,sizeof(int),1,fd);
	if(FileCounts!=NFiles)
	{
	  cout<<"File count mismatch for file "<<filename<<": expect "<<FileCounts<<", got "<<NFiles<<endl;
	  exit(1);
	}
	
	if(0==iFile)
	{
	  NumberOfHaloes=TotNgroups; 	  
	  Mvir.resize(TotNgroups);
	  Rvir.resize(TotNgroups);
	}
	
	if(IsGroupV3)
	fseek(fd, (sizeof(int)*2+sizeof(MyReal)*4)*Ngroups, SEEK_CUR);//skip len,offset,mass,pos
	else
	  fseek(fd, sizeof(int)*(2*Ngroups+3*Nsub), SEEK_CUR); 
	fseek(fd, sizeof(MyReal)*Ngroups*2*itype, SEEK_CUR);
	myfread(Mvir.data()+Nload, sizeof(MyReal), Ngroups, fd);
	myfread(Rvir.data()+Nload, sizeof(MyReal), Ngroups, fd);
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
  cout<<"Halo size "<<virtype<<"loaded at snapshot "<<SnapshotId<<", TotNgroups="<<TotNgroups<<" NFiles="<<NFiles<<endl;
  
}
}
#undef myfread

#endif
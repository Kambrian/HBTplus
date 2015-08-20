using namespace std;
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>

#include "simulation_io.h"
#include "../mymath.h"

#define myfread(a,b,c,d) fread_swap(a,b,c,d,ByteOrder)

void Snapshot_t::FormatSnapshotNumber(stringstream &ss)
{
  ss << setw(3) << setfill('0') << SnapshotIndex;
}
int Snapshot_t::GetSnapshotIndex()
{
  return SnapshotIndex;
}
void Snapshot_t::SetSnapshotIndex(int snapshot_index)
{
  SnapshotIndex=snapshot_index;
}
void Snapshot_t::CheckSnapshotIndexIsValid()
{
  if(SpecialConst::NullSnapshotId==SnapshotIndex)
  {
	cerr<<"Snapshot number not set!\n";
	exit(1);
  }
}
void Snapshot_t::GetSnapshotFileName(Parameter_t &param, int ifile, string &filename)
{
  FILE *fp;
  char buf[1024];
  sprintf(buf,"%s/snapdir_%03d/%s_%03d.%d",param.SnapshotPath,SnapshotIndex,param.SnapshotFileBase,SnapshotIndex,ifile);
  if(ifile==0)
	if(!file_exist(buf)) sprintf(buf,"%s/%s_%03d",param.SnapshotPath,param.SnapshotFileBase,SnapshotIndex); //try the other convention
  if(!file_exist(buf)) sprintf(buf,"%s/%s_%03d.%d",param.SnapshotPath,param.SnapshotFileBase,SnapshotIndex,ifile); //try the other convention
  if(!file_exist(buf)) sprintf(buf,"%s/%d/%s.%d",param.SnapshotPath,SnapshotIndex,param.SnapshotFileBase,ifile);//for BJL's RAMSES output
  if(!file_exist(buf))
  {
	cerr<<"Failed to find a snapshot file at "<<buf<<endl;
	exit(1);
  }
  filename=buf;
}
bool Snapshot_t::ReadFileHeader(FILE *fp, SnapshotHeader_t &header)
{
  //read the header part, assign header extensions, and check byteorder

  int headersize(SNAPSHOT_HEADER_SIZE),headersize_byteswap(SNAPSHOT_HEADER_SIZE);
  swap_Nbyte(&headersize_byteswap,1,sizeof(headersize_byteswap));

  bool ByteOrder;
  int dummy,dummy2;
  fread(&dummy,sizeof(dummy),1,fp);
  if(dummy==headersize)
	ByteOrder=false;
  else if(dummy==headersize_byteswap)
	ByteOrder=true;
  else
  {
	fprintf(stderr,"endianness check failed for header\n file format not expected:%d;%d,%d\n",dummy,headersize,headersize_byteswap);
	fflush(stderr);
	exit(1);
  }
  dummy=headersize;
  
  myfread(header.npart,sizeof(int),NUMBER_OF_PARTICLE_TYPES,fp);
  myfread(header.mass,sizeof(double),NUMBER_OF_PARTICLE_TYPES,fp);
  myfread(&header.time,sizeof(double),1,fp);
  myfread(&header.redshift,sizeof(double),1,fp);
  myfread(&header.flag_sfr,sizeof(int),1,fp);
  myfread(&header.flag_feedback,sizeof(int),1,fp);
  myfread(header.npartTotal,sizeof(int),NUMBER_OF_PARTICLE_TYPES,fp);
  myfread(&header.flag_cooling,sizeof(int),1,fp);
  myfread(&header.num_files,sizeof(int),1,fp);
  myfread(&header.BoxSize,sizeof(double),1,fp);
  myfread(&header.Omega0,sizeof(double),1,fp);
  myfread(&header.OmegaLambda,sizeof(double),1,fp);
  myfread(&header.HubbleParam,sizeof(double),1,fp);
  fseek(fp,headersize+sizeof(int),SEEK_SET);
  myfread(&dummy2,sizeof(dummy2),1,fp);
  if(dummy!=dummy2)
  {
	fprintf(stderr,"error!record brackets not match for header!\t%d,%d\n",dummy,dummy2);
	exit(1);
  } 
  
  /*extend and examine the header*/	  
  header.Hz=PhysicalConst::H0 * sqrt(header.Omega0 / (header.time * header.time * header.time) 
  + (1 - header.Omega0 - header.OmegaLambda) / (header.time * header.time)
  + header.OmegaLambda);//Hubble param for the current catalogue;
  
  return ByteOrder;
}

void Snapshot_t::LoadHeader(Parameter_t & param, int ifile=1)
{
  //read the header part, assign header extensions, and check byteorder
  CheckSnapshotIndexIsValid();

  FILE *fp; 
  string filename;
  GetSnapshotFileName(param, ifile, filename);
  myfopen(fp, filename, "r");
  
  ByteOrder=ReadFileHeader(fp, Header);
  
  if(Header.BoxSize!=param.BoxSize)
  {
	cerr<<"BoxSize not match input: read "<<Header.BoxSize<<"; expect "<<param.BoxSize<<endl;
	cerr<<"Maybe the length unit differ? Excpected unit: "<<param.LengthInMpch<<" Msol/h\n";
	exit(1);
  }
}

void Snapshot_t::Load(Parameter_t & param, bool load_id=true, bool load_position=true, bool load_velocity=true)
{
  LoadHeader(param);
  if(load_id)
	LoadId(param);
  if(load_position)
	LoadPosition(param);
  if(load_velocity)
	LoadVelocity(param);
}

void Snapshot_t::LoadId(Parameter_t &param)
{
  CheckSnapshotIndexIsValid();
//FIXME: finish this...........
  int dummy,dummy2;
  long pre_len,tail_len;
  HBTInt i,j;
  HBTInt Nload=0;
// 	void *PID=malloc(sizeof(IDatInt)*NP_DM);
	

  for(int iFile=0; iFile<NFILES; iFile++)
    {
	  FILE *fp; 
	  string filename;
	  GetSnapshotFileName(param, ifile, filename);
	  myfopen(fp, filename, "r");
	  
	  SnapshotHeader_t header;
	  bool ByteOrder=ReadFileHeader(fp, header);
	  
	  pre_len=header.npart[0]*sizeof(IDatReal)*3;
	  for(j=2,tail_len=0;j<6;j++)tail_len+=header.npart[j];
	  tail_len*=sizeof(IDatReal)*3;
      SKIP;
	  fseek(fp,pre_len,SEEK_CUR);//load only DM data;
	  if(flag_pos)
	  myfread(IDat.Pos+Nload,sizeof(IDatReal),3L*header.npart[1],fp);
	  else
	  fseek(fp,sizeof(IDatReal)*3*header.npart[1],SEEK_CUR);
	  fseek(fp,tail_len,SEEK_CUR);
      SKIP2;
	  CHECK(1);
	
      SKIP;
      fseek(fp,pre_len,SEEK_CUR);//load only DM data;
	  if(flag_vel)
		myfread(IDat.Vel+Nload,sizeof(IDatReal),3L*header.npart[1],fp);
	  else
		fseek(fp,sizeof(IDatReal)*3*header.npart[1],SEEK_CUR);
	  fseek(fp,tail_len,SEEK_CUR);
      SKIP2;
	  CHECK(2);
      
#ifndef NO_ID_RECORD //for Huiyuan's data.
	  pre_len=header.npart[0]*sizeof(IDatInt);
	  for(j=2,tail_len=0;j<6;j++)tail_len+=header.npart[j];
	  tail_len*=sizeof(IDatInt);
	  SKIP;
      fseek(fp,pre_len,SEEK_CUR);//load only DM data;
	  if(flag_id)
	  myfread(IDat.PID+Nload,sizeof(IDatInt),header.npart[1],fp);
	  else
	  fseek(fp,sizeof(IDatInt)*header.npart[1],SEEK_CUR);
	  fseek(fp,tail_len,SEEK_CUR);
      SKIP2;
	  CHECK(3);
#endif	  
	  
	  if(flag_scalar)
	  {
	  //gravity record
		SKIP;
		fseek(fp, sizeof(IDatReal)*3*NP_SIM, SEEK_CUR);
		SKIP2;
		CHECK(31);
	  //fifth force record
		SKIP;
		fseek(fp, sizeof(IDatReal)*3*NP_SIM, SEEK_CUR);
		SKIP2;
		CHECK(32);
	  //potential record
		pre_len=header.npart[0]*sizeof(IDatReal)*2;
		for(j=2,tail_len=0;j<6;j++)tail_len+=header.npart[j];
		tail_len*=sizeof(IDatReal)*2;
		SKIP;
		fseek(fp,pre_len,SEEK_CUR);//load only DM data;
		myfread(IDat.Pot+Nload,sizeof(IDatReal)*2,header.npart[1],fp);
		fseek(fp,tail_len,SEEK_CUR);
		SKIP2;
		CHECK(33);
	  }
		
	  if(feof(fp))
	  {
		fprintf(logfile,"error:End-of-File in %s\n",buf);
		fflush(logfile);exit(1);  
	  }
	
	  Nload+=header.npart[1];
//      printf("%d,%d,%ld\n",i, header.npart[1], Nload);fflush(stdout);
      fclose(fp);
  }
  if(NP_DM!=Nload)
  {
	  fprintf(logfile,"error: Number of loaded DM particles mismatch: %lld,%lld\n",(long long)NP_DM,(long long)Nload);
	  fflush(logfile);
	  exit(1);
  }

if(flag_id)
{
//now transfer to HBT's internal data  
  #ifdef HBTPID_RANKSTYLE
  struct io_ID2Ind *table;
  table=mymalloc(sizeof(struct io_ID2Ind)*NP_DM);
  for(i=0;i<NP_DM;i++)
  {
	  table[i].PID=IDat.PID[i];
	  table[i].PInd=i;
  }
  qsort(table,NP_DM,sizeof(struct io_ID2Ind),comp_PIDArr);
  Pdat.PID=mymalloc(sizeof(HBTInt)*NP_DM);
  for(i=0;i<NP_DM;i++)
	Pdat.PID[table[i].PInd]=i;  //now PID has been turned into particle ranks
	
  sprintf(buf,"%s/DM_PIDs_Sorted.dat",SUBCAT_DIR);
  if(!try_readfile(buf))//create the file if it does not exist
  {
	  IDatInt np;
	  myfopen(fp,buf,"w");
	  np=NP_DM;
	  fwrite(&np,sizeof(IDatInt),1,fp);
	  for(i=0;i<NP_DM;i++)
	  fwrite(&(table[i].PID),sizeof(IDatInt),1,fp);
	  fwrite(&np,sizeof(IDatInt),1,fp);
	  fclose(fp);
  }
	
  myfree(table);  
  myfree(IDat.PID);
  #else
#ifdef SAME_INTTYPE
  Pdat.PID=IDat.PID;
#else
  Pdat.PID=mymalloc(sizeof(HBTInt)*NP_DM);
  for(i=0;i<NP_DM;i++)
	Pdat.PID[i]=IDat.PID[i];
  myfree(IDat.PID);
#endif
  #endif
}
else
Pdat.PID=NULL;

if(flag_pos)
{
#ifdef SAME_REALTYPE
	  Pdat.Pos=IDat.Pos;
#else
	  Pdat.Pos=mymalloc(sizeof(HBTxyz)*NP_DM);
	  for(i=0;i<NP_DM;i++)
		for(j=0;j<3;j++)
			Pdat.Pos[i][j]=IDat.Pos[i][j];
	  myfree(IDat.Pos);
#endif
  	
  #ifdef CONVERT_LENGTH_MPC_KPC
  for(i=0;i<NP_DM;i++)
  for(j=0;j<3;j++)
  Pdat.Pos[i][j]*=1000.;
  #endif
  
  #ifdef PERIODIC_BDR
  for(i=0;i<NP_DM;i++)
  for(j=0;j<3;j++)
  Pdat.Pos[i][j]=position_modulus(Pdat.Pos[i][j]);
  #endif	
}
else
Pdat.Pos=NULL;

if(flag_vel)
{  
  #ifdef SAME_REALTYPE
	  Pdat.Vel=IDat.Vel;
#else
	  Pdat.Vel=mymalloc(sizeof(HBTxyz)*NP_DM);
	  for(i=0;i<NP_DM;i++)
		for(j=0;j<3;j++)
			Pdat.Vel[i][j]=IDat.Vel[i][j];
	  myfree(IDat.Vel);
#endif
}
else
Pdat.Vel=NULL;

if(flag_scalar)
{  
  Pdat.Vel=mymalloc(sizeof(HBTReal)*NP_DM);
  for(i=0;i<NP_DM;i++)
	Pdat.Scalar[i]=IDat.Pot[i][1];//pass scalar field to Pdat
  myfree(IDat.Pot);
}
else
  Pdat.Scalar=NULL;

}
void Snapshot_t::LoadPosition(Parameter_t &param)
{
  CheckSnapshotIndexIsValid();
}
void Snapshot_t::LoadVelocity(Parameter_t &param)
{
  CheckSnapshotIndexIsValid();
}

void Snapshot_t::Clear()
{
  delete [] ParticleId;
  delete [] ComovingPosition;
  delete [] PhysicalVelocity;
  NumberOfParticles=0;
}
HBTInt Snapshot_t::GetNumberOfParticles()
{
  return NumberOfParticles;
}
HBTInt Snapshot_t::GetParticleId(HBTInt index)
{
  return ParticleId[index];
}
HBTxyz& Snapshot_t::GetComovingPosition(HBTInt index)
{
  return ComovingPosition[index];
}
HBTxyz& Snapshot_t::GetPhysicalVelocity(HBTInt index)
{
  return PhysicalVelocity[index];
}

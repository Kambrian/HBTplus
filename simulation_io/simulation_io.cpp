using namespace std;
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "simulation_io.h"
#include "../mymath.h"

#define myfread(a,b,c,d) fread_swap(a,b,c,d,ByteOrder)

void Snapshot_t::FormatSnapshotNumber(stringstream &ss)
{
  ss << setw(3) << setfill('0') << Header.SnapshotIndex;
}
int Snapshot_t::GetSnapshotIndex()
{
  return Header.SnapshotIndex;
}
void Snapshot_t::SetSnapshotIndex(int snapshot_index)
{
  Header.SnapshotIndex=snapshot_index;
}
void Snapshot_t::CheckSnapshotIndexIsValid()
{
  if(SpecialConst::NullSnapshotId==Header.SnapshotIndex)
  {
	cerr<<"Snapshot number not initialized!\n";
	exit(1);
  }
}
void Snapshot_t::LoadHeader(char * snapshot_path, int iFile=0)
{
//read the header part, assign header extensions, and do several consistency check
//return ByteOrder
	int dummy,dummy2,n,ns,ByteOrder;
	
	FILE *fp;
	char buf[1024];

	sprintf(buf,"%s/snapdir_%03d/%s_%03d.%d",snapshot_path,(int)Nsnap,SNAPFILE_BASE,(int)Nsnap,ifile);
    if(1==NFILES)
	 if(!try_readfile(buf))	sprintf(buf,"%s/%s_%03d",snapshot_path,SNAPFILE_BASE,(int)Nsnap); //try the other convention
	
	if(!try_readfile(buf))	sprintf(buf,"%s/%s_%03d.%d",snapshot_path,SNAPFILE_BASE,(int)Nsnap,ifile); //try the other convention
	if(!try_readfile(buf))	sprintf(buf,"%s/%d/%s.%d",snapshot_path,(int)Nsnap,SNAPFILE_BASE,ifile);//for BJL's RAMSES output

	myfopen(fp,buf,"r");
	read_gadget_header(fp,h);
	fclose(fp);
	h->Nsnap=Nsnap;
	
	n=sizeof(IO_HEADER);
	ns=n;
	swap_Nbyte(&ns,1,sizeof(ns));
			
	fread(&dummy,sizeof(dummy),1,fp);
	
	if(dummy==n)
	 ByteOrder=0;
	else if(dummy==ns)
	 ByteOrder=1;
	else
	{
		fprintf(logfile,"endianness check failed for header\n file format not expected:%d;%d,%d\n",dummy,n,ns);
		fflush(logfile);
		exit(1);
	}
	
	dummy=n;

	myfread(h->npart,sizeof(int),6,fp);
	myfread(h->mass,sizeof(double),6,fp);
	myfread(&h->time,sizeof(double),1,fp);
	myfread(&h->redshift,sizeof(double),1,fp);
	myfread(&h->flag_sfr,sizeof(int),1,fp);
	myfread(&h->flag_feedback,sizeof(int),1,fp);
	myfread(h->npartTotal,sizeof(int),6,fp);
	myfread(&h->flag_cooling,sizeof(int),1,fp);
	myfread(&h->num_files,sizeof(int),1,fp);
	myfread(&h->BoxSize,sizeof(double),1,fp);
	myfread(&h->Omega0,sizeof(double),1,fp);
	myfread(&h->OmegaLambda,sizeof(double),1,fp);
	myfread(&h->HubbleParam,sizeof(double),1,fp);
	fseek(fp,n+sizeof(int),SEEK_SET);
	myfread(&dummy2,sizeof(dummy2),1,fp);
    if(dummy!=dummy2)
	{
		fprintf(logfile,"error!record brackets not match for header!\t%d,%d\n",dummy,dummy2);
		exit(1);
	} 
	
	/*extend and examine the header*/
	if(NFILES!=h->num_files)
	  {
		  fprintf(logfile,"error: number of snapfiles specified not the same as stored: %d,%d\n",
		  NFILES,h->num_files);
		  fflush(logfile);
//		  exit(1);
	  }
	  
	h->Hz=HUBBLE0 * sqrt(h->Omega0 / (h->time * h->time * h->time) 
			+ (1 - h->Omega0 - h->OmegaLambda) / (h->time * h->time)
			+ h->OmegaLambda);//Hubble param for the current catalogue;
	  
	#ifdef CONVERT_LENGTH_MPC_KPC
	  h->BoxSize*=1000.;
	#endif 
	
	return ByteOrder;
	

}
void Snapshot_t::Load(char * snapshot_path, bool load_id=true, bool load_position=true, bool load_velocity=true)
{
  if(load_id)
	LoadId();
  if(load_position)
	LoadPosition();
  if(load_velocity)
	LoadVelocity();
}
void Snapshot_t::LoadId()
{
  CheckSnapshotIndexIsValid();
}
void Snapshot_t::LoadPosition()
{
  CheckSnapshotIndexIsValid();
}
void Snapshot_t::LoadVelocity()
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

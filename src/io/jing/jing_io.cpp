#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstring>
#include <omp.h>

#include "../../snapshot.h"
#include "jing_io.h"
#include "fortread.h"

void JingReader_t::ReadParticleArray(float partarr[][3], long int np, int fileno)
{
  if(Header.ParticleDataXMajor)
    read_part_arr_xmajor_(partarr, &np, &fileno);
  else
    read_part_arr_imajor_(partarr, &np, &fileno);
}
void JingReader_t::ReadIdFileSingle(int ifile, vector< Particle_t >& Particles)
{
  long int nread=Header.Np/NumFilesId;
  string filename=GetFileName("id", ifile);
  vector <HBTInt> PId(nread);
  
  int fileno, filestat, flag_endian=NeedByteSwap;
  open_fortran_file_(filename.c_str(),&fileno,&flag_endian,&filestat);
  if(filestat) throw(runtime_error("failed to open file "+filename+", error no. "+to_string(filestat)+"\n"));
  read_fortran_record_HBTInt(PId.data(), &nread, &fileno);
  close_fortran_file_(&fileno);
  
  auto curr_particles=Particles.data()+nread*ifile;
  for(HBTInt i=0;i<nread;i++)
	curr_particles[i].Id=PId[i];
}
void JingReader_t::ReadId(vector <Particle_t> &Particles)
{
  if(HBTConfig.SnapshotHasIdBlock)
  {
    for(int ifile=0;ifile<NumFilesId;ifile++)
    #pragma omp task firstprivate(ifile) shared(Particles)
      ReadIdFileSingle(ifile, Particles);
  }
  else
  {
    #pragma omp task shared(Particles)
    for(HBTInt i=0;i<Particles.size();i++)
      Particles[i].Id=i;
  }
}
void JingReader_t::CheckIdRange(vector <Particle_t>&Particles)
{
  HBTInt i;
  #pragma omp parallel for
  for(i=0;i<Header.Np;i++)
    if(Particles[i].Id<1||Particles[i].Id>Header.Np)
      throw(runtime_error("id not in the range 1~"+to_string(Header.Np)+", for i="+to_string(i)+",pid="+to_string(Particles[i].Id)+", snap="+to_string(SnapshotId)+"\n"));
}
void JingReader_t::ReadPosFileSingle(int ifile, vector< Particle_t >& Particles)
{
  long int nread=Header.Np/NumFilesId;
  string filename=GetFileName("pos", ifile);
  auto curr_particles=Particles.data()+nread*ifile;
  
  int fileno,filestat, flag_endian=NeedByteSwap;
  open_fortran_file_(filename.c_str(),&fileno,&flag_endian,&filestat);
  if(filestat) throw(runtime_error("failed to open file "+filename+", error no. "+to_string(filestat)+"\n"));
	
  if(0==ifile)	skip_fortran_record_(&fileno);
	
  if(Header.FlagHasScale)//has scale
  {
    assert(sizeof(short)==2);
    vector <short> tmp(nread);
    for(int i=0;i<3;i++)
    {
      read_fortran_record2_(tmp.data(),&nread,&fileno);
      for(HBTInt j=0;j<nread;j++)
	      curr_particles[j].ComovingPosition[i]=tmp[j]*Header.xscale;
    }
  }
  else
  {
    typedef float PosVec_t[3];
    PosVec_t *Pos=new PosVec_t[nread];
    ReadParticleArray(Pos,nread,fileno);
    for(HBTInt j=0;j<nread;j++)
      copyHBTxyz(curr_particles[j].ComovingPosition,Pos[j]);
    delete [] Pos;
  }
  close_fortran_file_(&fileno);
}
void JingReader_t::ReadVelFileSingle(int ifile, vector< Particle_t >& Particles)
{
  long int nread=Header.Np/NumFilesId;
  string filename=GetFileName("vel", ifile);
  auto curr_particles=Particles.data()+nread*ifile;
  
  int fileno,filestat, flag_endian=NeedByteSwap;
  open_fortran_file_(filename.c_str(),&fileno,&flag_endian,&filestat);
  if(filestat) throw(runtime_error("failed to open file "+filename+", error no. "+to_string(filestat)+"\n"));
	
  if(0==ifile)	skip_fortran_record_(&fileno);
	
  if(Header.FlagHasScale)//has scale
  {
    assert(sizeof(short)==2);
    vector <short> tmp(nread);
    for(int i=0;i<3;i++)
    {
      read_fortran_record2_(tmp.data(),&nread,&fileno);
      for(HBTInt j=0;j<nread;j++)
	      curr_particles[j].PhysicalVelocity[i]=tmp[j]*Header.vscale;
    }
  }
  else
  {
    typedef float VelVec_t[3];
    VelVec_t *Vel=new VelVec_t[nread];
    ReadParticleArray(Vel,nread,fileno);
    for(HBTInt j=0;j<nread;j++)
      copyHBTxyz(curr_particles[j].PhysicalVelocity, Vel[j]);
    delete [] Vel;
  }
  close_fortran_file_(&fileno);
}
void JingReader_t::ReadPosition(vector <Particle_t> &Particles)
{
	for(int ifile=0;ifile<NumFilesPos;ifile++)
  	#pragma omp task firstprivate(ifile) shared(Particles)
	  ReadPosFileSingle(ifile, Particles);
}
void JingReader_t::ReadVelocity(vector <Particle_t> &Particles)
{
	for(int ifile=0;ifile<NumFilesVel;ifile++)
  	#pragma omp task firstprivate(ifile) shared(Particles)
	  ReadVelFileSingle(ifile, Particles);
}
string JingReader_t::GetFileName(const char * filetype, int iFile)
{
  char buf[1024];
  if(iFile==0)
  sprintf(buf,"%s/%s%s.%04d", HBTConfig.SnapshotPath.c_str(), filetype, HBTConfig.SnapshotFileBase.c_str(),SnapshotId);
  if(!file_exist(buf))  sprintf(buf,"%s/%s%s.%04d.%02d",HBTConfig.SnapshotPath.c_str(), filetype, HBTConfig.SnapshotFileBase.c_str(),SnapshotId, iFile+1);
  return string(buf);
}

void JingReader_t::ProbeFiles()
{ 
  int flag_endian=NeedByteSwap,filestat,fileno;
  int reset_fileno=1;
  alloc_file_unit_(&reset_fileno);//reset fortran fileno pool 
  
  string filename=GetFileName("pos");
  open_fortran_file_(filename.c_str(),&fileno,&flag_endian,&filestat);
  if(filestat)
  {
    cerr<<"Error opening file "<<filename<<",error no. "<<filestat<<endl;
    exit(1);
  }

  JingHeader_t header;
  read_part_header(&header.Np,&header.ips,&header.Redshift,&header.Omegat,&header.Lambdat,
					  &header.BoxSize,&header.xscale,&header.vscale,&fileno);
  auto L=header.BoxSize;
  if(header.ips!=SnapshotId)
  {
    NeedByteSwap=!NeedByteSwap;
    auto i=header.ips;
    swap_Nbyte(&i, 1, sizeof(i));
	swap_Nbyte(&L, 1, sizeof(L));
    if(i!=SnapshotId)
    {
      cerr<<"Error: fail to determine endianness for snapshot "<<SnapshotId<<", read ips="<<header.ips<<" or "<<i<<"(byteswapped)"<<endl;
      exit(1);
    }
  }
  if(L!=HBTConfig.BoxSize)
	cerr<<"Warning: boxsize does not match, maybe different units? "<<L<<","<<HBTConfig.BoxSize<<endl;

  close_fortran_file_(&fileno);

  NumFilesId=CountFiles("id");
  NumFilesPos=CountFiles("pos");
  NumFilesVel=CountFiles("vel");
}

int JingReader_t::CountFiles(const char *filetype)
{
  char filename[1024], pattern[1024], basefmt[1024], fmt[1024];
  const int ifile=1;
  sprintf(basefmt,"%s/%%s%s.%04d",HBTConfig.SnapshotPath.c_str(), HBTConfig.SnapshotFileBase.c_str(), SnapshotId);
  sprintf(basefmt, "%s/groups_%03d/subhalo_%%s_%03d",HBTConfig.HaloPath.c_str(),SnapshotId,SnapshotId);
  sprintf(fmt, "%s.%%02d", basefmt);
  sprintf(filename, fmt, filetype, ifile);
  int nfiles=1;
  if(file_exist(filename))
  {
	sprintf(pattern, basefmt, filetype);
	strcat(pattern, ".*");
	nfiles=count_pattern_files(pattern);
  }
  return nfiles;
}

void JingReader_t::ReadHeader(JingHeader_t& header, const char *filetype, int ifile)
{	
  short *tmp;
  int flag_endian=NeedByteSwap,filestat,fileno;
  int reset_fileno=1;
  alloc_file_unit_(&reset_fileno);//reset fortran fileno pool 
  
  string filename=GetFileName(filetype, ifile);
  open_fortran_file_(filename.c_str(),&fileno,&flag_endian,&filestat);
  if(filestat)
  {
    cerr<<"Error opening file "<<filename<<",error no. "<<filestat<<endl;
    exit(1);
  }
  read_part_header(&header.Np,&header.ips,&header.Redshift,&header.Omegat,&header.Lambdat,
		   &header.BoxSize,&header.xscale,&header.vscale,&fileno);
  close_fortran_file_(&fileno);
//   assert(header.BoxSize==HBTConfig.BoxSize);
  
  LoadExtraHeaderParams(header);
  header.mass[0]=0.;
  header.mass[1]=header.OmegaM0*3.*PhysicalConst::H0*PhysicalConst::H0/8./3.1415926/PhysicalConst::G*header.BoxSize*header.BoxSize*header.BoxSize/header.Np;//particle mass in units of 10^10Msun/h

  cout<<"z="<<header.Redshift<<endl;
    
  float Hratio; //(Hz/H0)
  float scale_reduced,scale0;//a,R0
  Hratio=sqrt(header.OmegaLambda0/header.Lambdat);
  header.Hz=PhysicalConst::H0*Hratio;
  scale_reduced=1./(1.+header.Redshift);
  header.ScaleFactor=scale_reduced;
  scale0=1+header.RedshiftIni;//scale_INI=1,scale_reduced_INI=1./(1.+z_ini),so scale0=scale_INI/scale_reduce_INI;
  header.vunit=100.*header.BoxSize*Hratio*scale_reduced*scale_reduced*scale0;   /*vunit=100*rLbox*R*(H*R)/(H0*R0)
  =L*H0*Hratio*R*R/R0 (H0=100 when length(L) in Mpc/h)
  *      =100*L*(H/H0)*a*a*R0
  * where a=R/R0;         */
}

void JingReader_t::LoadExtraHeaderParams(JingHeader_t &header)
{
  string filename=HBTConfig.ConfigFile;
  ifstream ifs;
  ifs.open(filename);
  if(!ifs.is_open())
  {
    cerr<<"Error opening parameter file "<<filename<<endl;
    exit(1);
  }
  string line;
  while(getline(ifs,line))
  {
    if(line.compare(0,13,"[ReaderExtra]")==0)//find the section
      break;
  }

  string name;
  HBTReal OmegaM0, OmegaL0, RedshiftIni;
  HBTInt SnapDivScale, ParticleDataXMajor;
#define ReadNameValue(x) ifs>>name>>x;assert(name==#x)
  ReadNameValue(OmegaM0);
  ReadNameValue(OmegaL0);
  ReadNameValue(RedshiftIni);
  ReadNameValue(SnapDivScale);
  ReadNameValue(ParticleDataXMajor);
#undef ReadNameValue
  ifs.close();
  
  header.OmegaM0=OmegaM0;
  header.OmegaLambda0=OmegaL0;
  header.RedshiftIni=RedshiftIni;
  header.SnapDivScale=SnapDivScale;
  header.FlagHasScale=(SnapshotId<=header.SnapDivScale);
  header.ParticleDataXMajor=ParticleDataXMajor;
}

void JingReader_t::LoadSnapshot(vector <Particle_t> &Particles, Cosmology_t &Cosmology)
{
  ReadHeader(Header);
  Cosmology.Set(Header.ScaleFactor, Header.OmegaM0, Header.OmegaLambda0);
#ifdef DM_ONLY
  Cosmology.ParticleMass=Header.mass[TypeDM];
#endif
  Particles.resize(Header.Np);
  
  int nthreads=NumFilesId+NumFilesPos+NumFilesVel;
  #pragma omp parallel num_threads(nthreads)
  #pragma omp single //nowait //creating tasks inside
  {
    ReadPosition(Particles);
    ReadVelocity(Particles);
    ReadId(Particles);
  }
// 	#pragma omp taskwait
	
  if(HBTConfig.SnapshotHasIdBlock)  CheckIdRange(Particles);
  
  #pragma omp parallel for
  for(HBTInt i=0;i<Header.Np;i++)
  {
    auto &p=Particles[i];
    for(int j=0;j<3;j++)
    {
      p.ComovingPosition[j]-=floor(p.ComovingPosition[j]);	//format coordinates to be in the range [0,1)
      p.ComovingPosition[j]*=Header.BoxSize;			//comoving coordinate in units of kpc/h
      p.PhysicalVelocity[j]*=Header.vunit;			//physical peculiar velocity in units of km/s
    }
  }
	
  cout<<" ( "<<nthreads<<" total files ) : "<<Particles.size()<<" particles loaded."<<endl;
}

namespace JingGroup
{
  bool check_b(float b)
  {
	return fabs(b-0.2)<0.01;
  }
  int ProbeGroupFileByteOrder(int snapshot_id)
  {
    int flag_endian=false, filestat, fileno;
    char buf[1024];
    sprintf(buf, "%s/fof.b20.%s.%04d",HBTConfig.HaloPath.c_str(), HBTConfig.SnapshotFileBase.c_str(), snapshot_id);
    open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
    if(filestat) throw(runtime_error("failed to open file "+string(buf)+", error no. "+to_string(filestat)+"\n"));
    float b;
    HBTInt Ngroups;
    read_group_header(&b,&Ngroups,&fileno);
    close_fortran_file_(&fileno);
    if(!check_b(b))
    {
      flag_endian=true;
      swap_Nbyte(&b, 1, sizeof(b));
      assert(check_b(b));
    }
    return flag_endian;
  }
  HBTInt LoadGroup(int snapshot_id, vector< Halo_t >& Halos)
  {
    long int nread;
	  
    int flag_endian=ProbeGroupFileByteOrder(snapshot_id), filestat, fileno;
    char buf[1024];
    sprintf(buf, "%s/fof.b20.%s.%04d",HBTConfig.HaloPath.c_str(), HBTConfig.SnapshotFileBase.c_str(), snapshot_id);
    open_fortran_file_(buf,&fileno,&flag_endian,&filestat);
    if(filestat) throw(runtime_error("failed to open file "+string(buf)+", error no. "+to_string(filestat)+"\n"));
    float b;
    HBTInt Ngroups;
    read_group_header(&b,&Ngroups,&fileno);
    Halos.resize(Ngroups);
    
    for(int i=0;i<3;i++)
      skip_fortran_record_(&fileno);
  //   vector <float> HaloCenX(Ngroups), HaloCenY(Ngroups), HaloCenZ(Ngroups);
  //   nread=Ngroups;
  //   read_fortran_record4_(HaloCenX.data(), &nread, &fileno);
  //   read_fortran_record4_(HaloCenY.data(), &nread, &fileno);
  //   read_fortran_record4_(HaloCenZ.data(), &nread, &fileno);
  //   for(HBTInt i=0;i<Ngroups;i++)
  //   {
  //     HaloCenX[i]*=HBTConfig.BoxSize;
  //     HaloCenY[i]*=HBTConfig.BoxSize;
  //     HaloCenZ[i]*=HBTConfig.BoxSize;
  //   }
  
    vector <HBTInt> Len(Ngroups), Offset(Ngroups);
    HBTInt Nids=0;
    for(HBTInt i=0;i<Ngroups;i++)
    {  
      nread=1;
      read_fortran_record_HBTInt(&Len[i],&nread,&fileno); 
      
      Offset[i]=Nids;
      Nids+=Len[i];
	  
      if(i>0&&Len[i]>Len[i-1]) 
	throw(runtime_error("Group size not ordered or wrong file format? group "+to_string(i)+", size "+to_string(Len[i])+", "+to_string(Len[i+1])+"\n"));
	  
      nread=Len[i];
      Halos[i].Particles.resize(Len[i]);
      read_fortran_record_HBTInt(Halos[i].Particles.data(),&nread,&fileno); 
    }
  
    close_fortran_file_(&fileno);
    
    cout<<"Snap="<<snapshot_id<<", Ngroups="<<Ngroups<<", Nids="<<Nids<<", b="<<b<<endl;
    
    if(HBTConfig.GroupLoadedIndex)
    #pragma omp parallel for
    for(HBTInt i=0;i<Ngroups;i++)
    {
	auto &particles=Halos[i].Particles;
	for(auto && p: particles)
	  p--;//change from [1,NP] to [0,NP-1] for index in C
    }
    
    return Nids;
  }
};


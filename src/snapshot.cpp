using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <chrono>
// #include <boost/concept_check.hpp>

#include "snapshot.h"
#include "mymath.h"
#include <mpi.h>

void create_Mpi_RemoteParticleType(MPI_Datatype& dtype, bool IdOnly)
{
  /*to create the struct data type for communication*/	
  RemoteParticle_t p;
  #define NumAttr 10
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint   offsets[NumAttr], origin,extent;
  
  MPI_Get_address(&p,&origin);
  MPI_Get_address((&p)+1,&extent);//to get the extent of s
  extent-=origin;
  
  int i=0;
  #define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
  RegisterAttr(Id, MPI_HBT_INT, 1)
  if(!IdOnly)
  {
  RegisterAttr(ComovingPosition, MPI_HBT_REAL, 3)
  RegisterAttr(PhysicalVelocity, MPI_HBT_REAL, 3)
  RegisterAttr(Mass, MPI_HBT_REAL, 1)
  }
  RegisterAttr(ProcessorId, MPI_INT, 1)
  RegisterAttr(Order, MPI_HBT_INT, 1)
  #undef RegisterAttr
  assert(i<=NumAttr);
  
  MPI_Type_create_struct(i,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef NumAttr
}
inline bool CompParticleProcessorId(const RemoteParticle_t &a, const RemoteParticle_t &b)
{
  return a.ProcessorId<b.ProcessorId;
}
inline bool CompParticleOrder(const RemoteParticle_t &a, const RemoteParticle_t &b)
{
  return a.Order<b.Order;
}

static void SortRemoteParticles(vector <RemoteParticle_t> &P)
{
  for(HBTInt i=0;i<P.size();i++)
  {
	auto &p=P[i];
	auto &j=p.Order;
	while(i!=j)
	  swap(p, P[j]);
  }
}

void ParticleExchanger_t::Exchange()
{
  //assumes all particles can be found on at least one processor
  auto it=ParticlesToProcess.begin();
  while(true)
  {
	bool alldone=ProcessParticle(it);
	if(alldone) break;
  }
  LocalParticles.assign(ParticlesToProcess.begin(), ParticlesToProcess.end());
  ParticlesToProcess.clear();
  
//   cout<<"remaining...\n";
  //send remaining particles
  while(ParticlesToSend.size())
	SendParticles();
  
  WaitChannel();
  
  RestoreParticles();
}

void ParticleExchanger_t::ReceiveParticles(bool blocking)
{
  int flag_ready=0;
  MPI_Status stat;
  if(blocking)
  {
	MPI_Probe(PrevRank, TagQuery, world.Communicator, &stat);
	flag_ready=1;
  }
  else
  {
	MPI_Iprobe(PrevRank, TagQuery, world.Communicator, &flag_ready, &stat);
// 	if(flag_ready) cout<<"Successful Iprobe!!\n";
  }

  if(flag_ready)
  {
// 	vector <RemoteParticle_t> recvbuffer(maxbuffersize);
	int buffersize;
	MPI_Get_count(&stat, MPI_RemoteParticleId_t, &buffersize);
	MPI_Recv(recvbuffer.data(), buffersize, MPI_RemoteParticleId_t, PrevRank, TagQuery, world.Communicator, MPI_STATUS_IGNORE);
	for(auto i=0;i<buffersize;i++)
	{
	  if(recvbuffer[i].Id==-1) iloop_debug++;
	  if(iloop_debug==nloop)
	  {
		if(i!=(buffersize-1))
		{
		  cout<<"bug: received message contains extra particles! "<<buffersize-1-i<<" ; total "<<buffersize<<endl;
		  for(int j=0;j<buffersize;j++){auto &p=recvbuffer[j]; cout<<"("<<p.Id<<","<<p.Order<<","<<p.ProcessorId<<")"<<endl;}
		  exit(1);
		}
	  }
	}
	ParticlesToProcess.insert(ParticlesToProcess.end(), recvbuffer.begin(), recvbuffer.begin()+buffersize);
  }
}

void ParticleExchanger_t::SendParticles()
{
	FlushChannel();
// 	vector <RemoteParticle_t> sendbuffer(maxbuffersize);
	int buffersize=min((HBTInt)maxbuffersize, (HBTInt)ParticlesToSend.size());
	if(ChannelIsClean&&buffersize)
	{
	  for(int i=0;i<buffersize;i++)
	  {
		sendbuffer[i]=ParticlesToSend.front();
		ParticlesToSend.pop_front();
	  }
	  for(int i=0;i<buffersize;i++)
	{
	  if(sendbuffer[i].Id==-1) iloop_sent++;
	  if(iloop_sent==nloop-1)
	  {
		if(i!=(buffersize-1))
		{
		  cout<<"bug: sent message contains extra particles! "<<buffersize-1-i<<" ; total "<<buffersize<<endl;
		  for(int j=0;j<buffersize;j++){auto &p=sendbuffer[j]; cout<<"("<<p.Id<<","<<p.Order<<","<<p.ProcessorId<<")"<<endl;}
		  exit(1);
		}
	  }
	}
	  MPI_Isend(sendbuffer.data(), buffersize, MPI_RemoteParticleId_t, NextRank, TagQuery, world.Communicator, &ReqSend);
	  ChannelIsClean=0;
	}
}
bool ParticleExchanger_t::ProcessParticle(ParticleStack_t::iterator &it)
{
	auto &p=*it;
	if(p.Id==EndParticleId)//end particle
	{
	  iloop++;
	  if(AllDone()) 
	  {
		it=ParticlesToProcess.erase(it);
		while(it!=ParticlesToProcess.end())
		{
		  cout<<"("<<it->Id<<","<<it->ProcessorId<<","<<it->Order<<": "<<world.rank()<<")";
		  ++it;
		}
		assert(it==ParticlesToProcess.end());//all particles should have been processed by now, unless some particles cannot be located on any processor!
		return true;
	  }
	  ParticlesToSend.push_back(p);
	  it=ParticlesToProcess.erase(it);
	}
	else
	{
	  HBTInt ind=snap.GetIndex(p);
	  if(ind!=SpecialConst::NullParticleId)
	  {
		p=snap.Particles[ind];
		++it;
	  }
	  else
	  {
		ParticlesToSend.push_back(p);
		it=ParticlesToProcess.erase(it);
	  }
	}
	bool is_starving=(it==ParticlesToProcess.end());
	if(is_starving)
	{
	  SendParticles();
	  int flag_begin=0;
	  if(it==ParticlesToProcess.begin())
		flag_begin=1;
	  else
		--it;
	  
	  ReceiveParticles(1);//blocking receive
	  
	  if(flag_begin)
		it=ParticlesToProcess.begin();
	  else
		++it;//to ensure iterator validity
	}
	else
	{
	  if(ParticlesToSend.size()>=maxbuffersize)
		SendParticles();
	  ReceiveParticles(0);//non-blocking
	}
// 	FlushChannel();//alternatively, can flush only inside send
	return false;
}

void ParticleExchanger_t::RestoreParticles()
{//move LocalParticles to particles on the original processor
  typedef vector <RemoteParticle_t>::iterator InputIterator_t;
  typedef vector <RemoteParticle_t>::iterator OutputIterator_t;
  sort(LocalParticles.begin(), LocalParticles.end(), CompParticleProcessorId);  //not necessary, should be in order already.
  vector <HBTInt> SendSizes(world.size(),0), SendDisps(world.size()), RecvSizes(world.size()), RecvDisps(world.size());
  for(auto &&p: LocalParticles)
	SendSizes[p.ProcessorId]++;
  CompileOffsets(SendSizes, SendDisps);
  MPI_Alltoall(SendSizes.data(), 1, MPI_HBT_INT, RecvSizes.data(), 1, MPI_HBT_INT, world.Communicator);
  HBTInt nrecv=CompileOffsets(RecvSizes, RecvDisps);
  vector <RemoteParticle_t> particles(nrecv);
  vector <InputIterator_t > SendBegin(world.size());
  vector <OutputIterator_t> RecvBegin(world.size());
  for(int i=0;i<world.size();i++)
  {
	SendBegin[i]=LocalParticles.begin()+SendDisps[i];
	RecvBegin[i]=particles.begin()+RecvDisps[i];
  }
  MyAllToAll<RemoteParticle_t, InputIterator_t, OutputIterator_t>(world, SendBegin, SendSizes, RecvBegin, MPI_RemoteParticle_t);
  for(int i=0;i<world.size();i++)
	RecvBegin[i]=particles.begin()+RecvDisps[i];
  for(HBTInt rank=0;rank<world.size();rank++)
  {
	for(int i=0;i<RecvSizes[rank];i++)
	  RecvBegin[rank][i].ProcessorId=rank;
  }
  SortRemoteParticles(particles);
  LocalParticles.swap(particles);
}

void SnapshotHeader_t::create_MPI_type(MPI_Datatype& dtype)
{
  /*to create the struct data type for communication*/	
  SnapshotHeader_t &p=*this;
  #define NumAttr 13
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint   offsets[NumAttr], origin,extent;
  
  MPI_Get_address(&p,&origin);
  MPI_Get_address((&p)+1,&extent);//to get the extent of s
  extent-=origin;
  
  int i=0;
  #define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
  RegisterAttr(npart, MPI_INT, NUMBER_OF_PARTICLE_TYPES)
  RegisterAttr(mass, MPI_DOUBLE, NUMBER_OF_PARTICLE_TYPES)
  RegisterAttr(ScaleFactor, MPI_DOUBLE, 1)
  RegisterAttr(redshift, MPI_DOUBLE, 1)
  RegisterAttr(flag_sfr, MPI_INT, 1)
  RegisterAttr(flag_feedback, MPI_INT, 1)
  RegisterAttr(npartTotal, MPI_UNSIGNED, 1)
  RegisterAttr(flag_cooling, MPI_INT, 1)
  RegisterAttr(num_files, MPI_INT, 1)
  RegisterAttr(BoxSize, MPI_DOUBLE, 1)
  RegisterAttr(OmegaM0, MPI_DOUBLE, 1)
  RegisterAttr(OmegaLambda0, MPI_DOUBLE, 1)
  RegisterAttr(HubbleParam, MPI_DOUBLE, 1)
  #undef RegisterAttr
  assert(i==NumAttr);
  
  MPI_Type_create_struct(NumAttr,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef NumAttr
}
ostream& operator << (ostream& o, Particle_t &p)
{
   o << "[" << p.Id << ", " <<p.Mass<<", " << p.ComovingPosition << ", " << p.PhysicalVelocity << "]";
   return o;
};
/*
void ParticleSnapshot_t::FillParticleHash()
{
  cout<<"Filling Hash Table...\n";
  auto begin = chrono::high_resolution_clock::now();
  ParticleHash.rehash(NumberOfParticles);
  ParticleHash.reserve(NumberOfParticles);
  for(HBTInt i=0;i<NumberOfParticles;i++)
	ParticleHash[ParticleId[i]]=i;
  auto end = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::milliseconds>(end - begin);
  cout << "inserts: " << elapsed.count() << endl;
//   cout<<ParticleHash.bucket_count()<<" buckets used; load factor "<<ParticleHash.load_factor()<<endl;
}
void ParticleSnapshot_t::ClearParticleHash()
{
  ParticleHash.clear();
}
*/
class ParticleKeyList_t: public KeyList_t <HBTInt, HBTInt>
{
  typedef HBTInt Index_t;
  typedef HBTInt Key_t;
  const ParticleSnapshot_t &Snap;
public:
  ParticleKeyList_t(ParticleSnapshot_t &snap): Snap(snap) {};
  Index_t size() const
  {
	return Snap.size();
  }
  Key_t GetKey(Index_t i) const
  {
	return Snap.GetId(i);
  }
  Index_t GetIndex(Index_t i) const
  {
	return i;
  }
};

void ParticleSnapshot_t::FillParticleHash()
{
  ParticleKeyList_t Ids(*this); 
  ParticleHash->Fill(Ids);
}
void ParticleSnapshot_t::ClearParticleHash()
{
  ParticleHash->Clear();
}

void ParticleSnapshot_t::AveragePosition(HBTxyz& CoM, const HBTInt Particles[], HBTInt NumPart) const
/*mass weighted average position*/
{
	HBTInt i,j;
	double sx[3],origin[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoM, GetComovingPosition(Particles[0]));
	  return;
	}
	
	sx[0]=sx[1]=sx[2]=0.;
	msum=0.;
	if(HBTConfig.PeriodicBoundaryOn)
	  for(j=0;j<3;j++)
		origin[j]=GetComovingPosition(Particles[0])[j];
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=GetMass(Particles[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
	  if(HBTConfig.PeriodicBoundaryOn)
		  sx[j]+=NEAREST(GetComovingPosition(Particles[i])[j]-origin[j])*m;
	  else
		  sx[j]+=GetComovingPosition(Particles[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	{
		sx[j]/=msum;
		if(HBTConfig.PeriodicBoundaryOn) sx[j]+=origin[j];
		CoM[j]=sx[j];
	}
}
void ParticleSnapshot_t::AverageVelocity(HBTxyz& CoV, const HBTInt Particles[], HBTInt NumPart) const
/*mass weighted average velocity*/
{
	HBTInt i,j;
	double sv[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoV, GetPhysicalVelocity(Particles[0]));
	  return;
	}
	
	sv[0]=sv[1]=sv[2]=0.;
	msum=0.;
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=GetMass(Particles[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
		sv[j]+=GetPhysicalVelocity(Particles[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	  CoV[j]=sv[j]/msum;
}

void Particle_t::create_MPI_type(MPI_Datatype &MPI_HBTParticle_t)
{
/*to create the struct data type for communication*/	
Particle_t &p=*this;
#define NumAttr 4
MPI_Datatype oldtypes[NumAttr];
MPI_Aint   offsets[NumAttr], origin,extent;
int blockcounts[NumAttr];

MPI_Get_address(&p,&origin);
MPI_Get_address(&p.Id,offsets);
MPI_Get_address(p.ComovingPosition.data(),offsets+1);//caution: this might be implementation dependent??
MPI_Get_address(p.PhysicalVelocity.data(),offsets+2);
MPI_Get_address(&p.Mass,offsets+3);
MPI_Get_address((&p)+1,&extent);//to get the extent of s

for(int i=0;i<NumAttr;i++)
  offsets[i]-=origin;

oldtypes[0] = MPI_HBT_INT;
blockcounts[0] = 1;

oldtypes[1] = MPI_HBT_REAL;
blockcounts[1] = 3;

oldtypes[2] = MPI_HBT_REAL;
blockcounts[2] = 3;

oldtypes[3] = MPI_HBT_REAL;
blockcounts[3] = 1;

extent-=origin;

assert(offsets[2]-offsets[1]==sizeof(HBTReal)*3);//to make sure HBTxyz is stored locally.
/*
for(int i=0;i<NumAttr;i++)
  cout<<offsets[i]<<',';
cout<<endl<<extent<<endl;
*/
MPI_Type_create_struct(NumAttr,blockcounts,offsets,oldtypes, &MPI_HBTParticle_t);//some padding is added automatically by MPI as well
MPI_Type_create_resized(MPI_HBTParticle_t,(MPI_Aint)0, extent, &MPI_HBTParticle_t);
MPI_Type_commit(&MPI_HBTParticle_t);
//~ MPI_Type_get_extent(*pMPIshearprof,&origin,&extent);
//~ printf("%d\n",extent);
#undef NumAttr
}

void AveragePosition(HBTxyz& CoM, const Particle_t Particles[], HBTInt NumPart)
/*mass weighted average position*/
{
	HBTInt i,j;
	double sx[3],origin[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoM, Particles[0].ComovingPosition);
	  return;
	}
	
	sx[0]=sx[1]=sx[2]=0.;
	msum=0.;
	if(HBTConfig.PeriodicBoundaryOn)
	  for(j=0;j<3;j++)
		origin[j]=Particles[0].ComovingPosition[j];
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=Particles[i].Mass;
	  msum+=m;
	  for(j=0;j<3;j++)
	  if(HBTConfig.PeriodicBoundaryOn)
		  sx[j]+=NEAREST(Particles[i].ComovingPosition[j]-origin[j])*m;
	  else
		  sx[j]+=Particles[i].ComovingPosition[j]*m;
	}
	
	for(j=0;j<3;j++)
	{
		sx[j]/=msum;
		if(HBTConfig.PeriodicBoundaryOn) sx[j]+=origin[j];
		CoM[j]=sx[j];
	}
}
void AverageVelocity(HBTxyz& CoV, const Particle_t Particles[], HBTInt NumPart)
/*mass weighted average velocity*/
{
	HBTInt i,j;
	double sv[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoV, Particles[0].PhysicalVelocity);
	  return;
	}
	
	sv[0]=sv[1]=sv[2]=0.;
	msum=0.;
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=Particles[i].Mass;
	  msum+=m;
	  for(j=0;j<3;j++)
		sv[j]+=Particles[i].PhysicalVelocity[j]*m;
	}
	
	for(j=0;j<3;j++)
	  CoV[j]=sv[j]/msum;
}

void Snapshot_t::SphericalOverdensitySize(float& Mvir, float& Rvir, HBTReal VirialFactor, const vector< HBTReal >& RSorted, HBTReal ParticleMass) const
/*
 * find SphericalOverdensitySize from a given list of sorted particle distances.
 * all distances comoving.
 * 
 * Brute-force method to scan from most-distant particle.
 * 
 * currently only support constant particle mass.
 */
{  
  HBTInt i, np=RSorted.size();
  HBTReal RhoVirial=VirialFactor*Hz*Hz/2.0/PhysicalConst::G/ParticleMass*ScaleFactor*ScaleFactor*ScaleFactor;
  for(i=np;i>0;i--)
  {
	HBTReal r=RSorted[i-1];
	if(i>RhoVirial*r*r*r) break;
  }
  Mvir=i*ParticleMass;
  Rvir=pow(i/RhoVirial, 1.0/3);//comoving
}

void Snapshot_t::SphericalOverdensitySize2(float& Mvir, float& Rvir, HBTReal VirialFactor, const vector< HBTReal >& RSorted, HBTReal ParticleMass) const
/*
 * find SphericalOverdensitySize from a given list of sorted particle distances.
 * all distances comoving.
 * 
 * Iterative method to guess the radius.
 * 
 * currently only support constant particle mass.
 */
{
  HBTReal tol=1e-5;
  HBTInt i,ndiv, np=RSorted.size();
  HBTReal rvir,rdiv;
  
  ndiv=np;//guess mass
  rdiv=RSorted[ndiv-1];
  rvir=pow(2.0*PhysicalConst::G*ndiv*ParticleMass/VirialFactor/Hz/Hz,1.0/3)/ScaleFactor;//guess radius
  while(true)
  {
	if(rdiv>rvir)//reduce mass guess
	{
	  for(i=ndiv-1;i>=0;i--)
	  {
		if(RSorted[i]<rvir) break;
	  }
	  ndiv=i+1;
	}
	else if(rdiv<rvir) //increase mass guess
	{
	  for(i=ndiv;i<np;i++)
	  {
		if(RSorted[i]>=rvir) break;
	  }
	  ndiv=i;
	}
	
	rdiv=rvir;
	rvir=pow(2.0*PhysicalConst::G*ndiv*ParticleMass/VirialFactor/Hz/Hz,1.0/3)/ScaleFactor;//recalibrate radius
	
	if(0==ndiv||np==ndiv||fabs((rvir-rdiv)/rvir)<tol) break; //converged
  }
  
  Rvir=rvir;
  Mvir=ndiv*ParticleMass;
}

void Snapshot_t::HaloVirialFactors(HBTReal &virialF_tophat, HBTReal &virialF_b200, HBTReal &virialF_c200) const
{
	HBTReal Hratio,x,OmegaZ;
	Hratio=Hz/PhysicalConst::H0;
	OmegaZ=OmegaM0/(ScaleFactor*ScaleFactor*ScaleFactor)/Hratio/Hratio;
	x=OmegaZ-1;
	virialF_tophat=18.0*3.1416*3.1416+82.0*x-39.0*x*x;//<Rho_vir>/Rho_cri
	virialF_c200=200.;
	virialF_b200=200.*OmegaZ;//virialF w.r.t contemporary critical density 
}
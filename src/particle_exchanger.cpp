#include "snapshot.h"
#include "particle_exchanger.h"

void create_Mpi_RemoteParticleType(MPI_Datatype& dtype)
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
  RegisterAttr(ComovingPosition, MPI_HBT_REAL, 3)
  RegisterAttr(PhysicalVelocity, MPI_HBT_REAL, 3)
  RegisterAttr(Mass, MPI_HBT_REAL, 1)
  RegisterAttr(ProcessorId, MPI_INT, 1)
  RegisterAttr(Order, MPI_HBT_INT, 1)
  #undef RegisterAttr
  assert(i<=NumAttr);
  
  MPI_Type_create_struct(i,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef NumAttr
}

void create_Mpi_RemoteParticleIdType(MPI_Datatype& dtype)
{
  /*to create the struct data type for communication*/	
  RemoteParticleId_t p;
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
inline bool CompParticleId(const RemoteParticleId_t &a, const RemoteParticleId_t &b)
{
  return a.Id<b.Id;
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

void ParticleExchanger_t::BcastParticles(HBTInt& ParticleCount)
{
  int root=CurrSendingRank;
  MPI_Bcast(&ParticleCount, 1, MPI_HBT_INT, root, world.Communicator);
  if(ParticleCount)
  {
	//determine loops
	const int chunksize=1024*1024;
	HBTInt  Nloop=ceil(1.*ParticleCount/chunksize);
	int buffersize=ParticleCount/Nloop+1, nremainder=ParticleCount%Nloop;
	//transmit
	vector <RemoteParticleId_t> buffer(buffersize);
	for(HBTInt iloop=0;iloop<Nloop;iloop++)
	{
	  if(iloop==nremainder)//switch sendcount from n+1 to n
	  {
		buffersize--;
		buffer.resize(buffersize);
	  }
	  if(world.rank()==root)//pack
	  {
		for(auto it_buff=buffer.begin();it_buff!=buffer.end();++it_buff)
		{
		  *it_buff=ParticlesToSend.back();
		  ParticlesToSend.pop_back();
		}
	  }
	  MPI_Bcast(buffer.data(), buffersize, MPI_RemoteParticleId_t, root, world.Communicator);
	  for(auto it_buff=buffer.begin();it_buff!=buffer.end();++it_buff)//unpack
		ParticlesToProcess.push_back(*it_buff);
	}
  }
  if(world.rank()==root)
  {
	SendStackSize-=ParticleCount;
	if(SendStackSize==0) CurrSendingRank++;
  }
  MPI_Bcast(&CurrSendingRank, 1, MPI_HBT_INT, root, world.Communicator);
}

bool ParticleExchanger_t::GatherParticles(HBTInt capacity)
{
  ParticlesToProcess.reserve(capacity);
  HBTInt nsend;
  while(capacity)
  {
	if(world.rank()==CurrSendingRank)
	  nsend=min(SendStackSize, capacity);
	BcastParticles(nsend);
	capacity-=nsend;
	if(CurrSendingRank==world.size()) return true;
  }
  
  return false;
}
void ParticleExchanger_t::QueryParticles()
{
  sort(ParticlesToProcess.begin(), ParticlesToProcess.end(), CompParticleId);

  snap.GetIndices(ParticlesToProcess);

  for(auto &&p: ParticlesToProcess)
  {
	if(p.Id!=SpecialConst::NullParticleId)
	  LocalParticles.emplace_back(snap.Particles[p.Id], p.ProcessorId, p.Order);
  }
 
 ParticlesToProcess.clear();
}

void ParticleExchanger_t::Exchange()
{
  HBTInt capacity=ceil(1.*snap.NumberOfParticlesOnAllNodes/world.size());//decrease this if out of memory; increase this to increase efficiency
  LocalParticles.reserve(capacity);
  while(true)
  {
	bool flag_end=GatherParticles(capacity);
	QueryParticles();//Walk through, append to LocalParticles, clear ParticlesToProcess
	if(flag_end) break;
  }

  RestoreParticles();
}

void ParticleExchanger_t::RestoreParticles()
{//move LocalParticles to particles on the original processor
  typedef vector <RemoteParticle_t>::iterator InputIterator_t;
  typedef vector <RemoteParticle_t>::iterator OutputIterator_t;
  
  sort(LocalParticles.begin(), LocalParticles.end(), CompParticleProcessorId);  
  
  vector <HBTInt> SendSizes(world.size(),0), SendDisps(world.size()), RecvSizes(world.size()), RecvDisps(world.size());
  for(auto &&p: LocalParticles)
	SendSizes[p.ProcessorId]++;
  CompileOffsets(SendSizes, SendDisps);
  MPI_Alltoall(SendSizes.data(), 1, MPI_HBT_INT, RecvSizes.data(), 1, MPI_HBT_INT, world.Communicator);
  HBTInt nrecv=CompileOffsets(RecvSizes, RecvDisps);
  assert(nrecv==SendStackSize0);
  
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
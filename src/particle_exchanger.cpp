#include "snapshot.h"
#include "particle_exchanger.h"

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
  
  OpenChannel();
  while(true)
  {
	bool alldone=ProcessParticle(it);
	if(alldone) break;
  }
  CloseChannel();
 
  LocalParticles.assign(ParticlesToProcess.begin(), ParticlesToProcess.end());
  ParticlesToProcess.clear();
  
//   cout<<"remaining...\n";
  //send remaining particles
  while(SendStackSize)
	SendParticles();
  
  RestoreParticles();
}

void ParticleExchanger_t::ReceiveParticles(bool blocking)
{
  int MessageArrived=0;
  MPI_Status bufferstat;
  if(blocking)
  {
	MPI_Wait(&ReqRecv, &bufferstat);
	MessageArrived=1;
  }
  else
	MPI_Test(&ReqRecv, &MessageArrived, &bufferstat);
  
  if(MessageArrived)
  {
	int buffersize;
	MPI_Get_count(&bufferstat, MPI_RemoteParticleId_t, &buffersize);
	ParticlesToProcess.insert(ParticlesToProcess.end(), RecvBuffer.begin(), RecvBuffer.begin()+buffersize);
	OpenChannel();
  }
}

void ParticleExchanger_t::SendParticles()
{
	int buffersize=min(maxbuffersize, SendStackSize);
	if(buffersize)
	{
	  for(int i=0;i<buffersize;i++)
	  {
		SendBuffer[i]=ParticlesToSend.front();
		ParticlesToSend.pop_front();
	  }
	  SendStackSize-=buffersize;
	  MPI_Send(SendBuffer.data(), buffersize, MPI_RemoteParticleId_t, NextRank, TagQuery, world.Communicator);
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
		assert(it==ParticlesToProcess.end());//all particles should have been processed by now, unless some particles cannot be located on any processor!
		return true;
	  }
	  ParticlesToSend.push_back(p);
	  SendStackSize++;
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
		SendStackSize++;
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
	  if(SendStackSize>=maxbuffersize)
		SendParticles();
	  ReceiveParticles(0);//non-blocking
	}
	return false;
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
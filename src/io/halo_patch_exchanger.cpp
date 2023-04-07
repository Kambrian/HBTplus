
#include <assert.h>
#include <list>

#include "../mymath.h"
#include "../halo.h"
#include "halo_patch_exchanger.h"

namespace HaloPatchExchanger
{

struct HaloInfo_t
{
  HBTInt id;
  HBTReal m;
  HBTxyz x;
  int order;
};

static void create_MPI_HaloInfo_t(MPI_Datatype &dtype)
{
  HaloInfo_t p;
  #define NumAttr 13
  MPI_Datatype oldtypes[NumAttr];
  int blockcounts[NumAttr];
  MPI_Aint   offsets[NumAttr], origin,extent;

  MPI_Get_address(&p,&origin);
  MPI_Get_address((&p)+1,&extent);//to get the extent of s
  extent-=origin;

  int i=0;
  #define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
  RegisterAttr(id, MPI_HBT_INT, 1)
  RegisterAttr(m, MPI_HBT_REAL, 1)
  RegisterAttr(x[0], MPI_HBT_REAL, 3)
  RegisterAttr(order, MPI_INT, 1)
  #undef RegisterAttr
  assert(i<=NumAttr);

  MPI_Type_create_struct(i,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef NumAttr
}

inline bool CompHaloInfo_Id(const HaloInfo_t &a, const HaloInfo_t &b)
{
  return a.id<b.id;
}
inline bool CompHaloInfo_Order(const HaloInfo_t &a, const HaloInfo_t &b)
{
  return a.order<b.order;
}
inline bool CompHaloId(const Halo_t &a, const Halo_t &b)
{
  return a.HaloId<b.HaloId;
}
double ReduceHaloPosition(vector <HaloInfo_t>::iterator it_begin, vector <HaloInfo_t>::iterator it_end, HBTxyz &x)
{
  HBTInt j;
  double sx[3],origin[3],msum;

  if(it_begin==it_end) return 0.;
  if(it_begin+1==it_end)
  {
    copyHBTxyz(x, it_begin->x);
    return it_begin->m;
  }

  sx[0]=sx[1]=sx[2]=0.;
  msum=0.;
  if(HBTConfig.PeriodicBoundaryOn)
    for(j=0;j<3;j++)
	  origin[j]=it_begin->x[j];

  for(auto it=it_begin;it!=it_end;++it)
  {
    HBTReal m=it->m;
    msum+=m;
    for(j=0;j<3;j++)
    if(HBTConfig.PeriodicBoundaryOn)
	    sx[j]+=NEAREST(it->x[j]-origin[j])*m;
    else
	    sx[j]+=it->x[j]*m;
  }

  for(j=0;j<3;j++)
  {
	  sx[j]/=msum;
	  if(HBTConfig.PeriodicBoundaryOn)
	  {
	    sx[j]+=origin[j];
	    x[j]=position_modulus(sx[j], HBTConfig.BoxSize);
	  }
	  else
	    x[j]=sx[j];
  }
  return msum;
}
void ReduceHaloRank(vector <HaloInfo_t>::iterator it_begin, vector <HaloInfo_t>::iterator it_end, HBTxyz &step, vector <int> &dims)
{
  HBTxyz x;
  ReduceHaloPosition(it_begin, it_end,x);
  int rank=AssignCell(x, step, dims);
  for(auto it=it_begin;it!=it_end;++it)
    it->id=rank; //store destination rank in id.
}
static void DecideTargetProcessor(MpiWorker_t& world, vector< Halo_t >& Halos, vector <IdRank_t> &TargetRank)
{
  int this_rank=world.rank();
  for(auto &&h: Halos)
    h.Mass=AveragePosition(h.ComovingAveragePosition, h.Particles.data(), h.Particles.size());

  vector <HaloInfo_t> HaloInfoSend(Halos.size()), HaloInfoRecv;
  for(HBTInt i=0;i<Halos.size();i++)
  {
    HaloInfoSend[i].id=Halos[i].HaloId;
    HaloInfoSend[i].m=Halos[i].Mass;
    HaloInfoSend[i].x=Halos[i].ComovingAveragePosition;
//     HaloInfoSend[i].rank=this_rank;
  }
  HBTInt MaxHaloId=0;
  if(Halos.size()) MaxHaloId=Halos.back().HaloId;
  MPI_Allreduce(MPI_IN_PLACE, &MaxHaloId, 1, MPI_HBT_INT, MPI_MAX, world.Communicator);
  HBTInt ndiv=(++MaxHaloId)/world.size();
  if(MaxHaloId%world.size()) ndiv++;
  vector <int> SendSizes(world.size(),0), SendOffsets(world.size()), RecvSizes(world.size()), RecvOffsets(world.size());
  for(HBTInt i=0;i<Halos.size();i++)
  {
    int idiv=Halos[i].HaloId/ndiv;
    SendSizes[idiv]++;
  }
  CompileOffsets(SendSizes, SendOffsets);
  MPI_Alltoall(SendSizes.data(), 1, MPI_INT, RecvSizes.data(), 1, MPI_INT, world.Communicator);
  int nhalo_recv=CompileOffsets(RecvSizes, RecvOffsets);
  HaloInfoRecv.resize(nhalo_recv);
  MPI_Datatype MPI_HaloInfo_t;
  create_MPI_HaloInfo_t(MPI_HaloInfo_t);
  MPI_Alltoallv(HaloInfoSend.data(), SendSizes.data(), SendOffsets.data(), MPI_HaloInfo_t, HaloInfoRecv.data(), RecvSizes.data(), RecvOffsets.data(), MPI_HaloInfo_t, world.Communicator);
  for(int i=0;i<nhalo_recv;i++)
    HaloInfoRecv[i].order=i;
  sort(HaloInfoRecv.begin(), HaloInfoRecv.end(), CompHaloInfo_Id);
  list <int> haloid_offsets;
  HBTInt curr_id=-1;
  for(int i=0;i<nhalo_recv;i++)
  {
    if(curr_id!=HaloInfoRecv[i].id)
    {
      haloid_offsets.push_back(i);
      curr_id=HaloInfoRecv[i].id;
    }
  }
  haloid_offsets.push_back(nhalo_recv);
  //combine coordinates and determine target
  auto dims=ClosestFactors(world.size(), 3);
  HBTxyz step;
  for(int i=0;i<3;i++)
	step[i]=HBTConfig.BoxSize/dims[i];
  auto it_end=haloid_offsets.end(); --it_end;
  for(auto it=haloid_offsets.begin();it!=it_end;it++)
  {
    auto it_next=it;
    ++it_next;
    ReduceHaloRank(HaloInfoRecv.begin()+*it, HaloInfoRecv.begin()+*it_next, step, dims);
  }
  sort(HaloInfoRecv.begin(), HaloInfoRecv.end(), CompHaloInfo_Order);
  //send back
  MPI_Alltoallv(HaloInfoRecv.data(), RecvSizes.data(), RecvOffsets.data(), MPI_HaloInfo_t, HaloInfoSend.data(), SendSizes.data(), SendOffsets.data(), MPI_HaloInfo_t, world.Communicator);
  MPI_Type_free(&MPI_HaloInfo_t);

  TargetRank.resize(Halos.size());
  for(HBTInt i=0; i<TargetRank.size();i++)
  {
    TargetRank[i].Id=i;
    TargetRank[i].Rank=HaloInfoSend[i].id;
  }

}

void MergeHalos(vector< Halo_t >& Halos)
{
  if(Halos.empty()) return;
  sort(Halos.begin(), Halos.end(), CompHaloId);
  auto it1=Halos.begin();
  for(auto it2=it1+1;it2!=Halos.end();++it2)
  {
    if(it2->HaloId==it1->HaloId)
    {
      it1->Particles.insert(it1->Particles.end(), it2->Particles.begin(), it2->Particles.end());
    }
    else
    {
      ++it1;
      if(it2!=it1)
	*it1=move(*it2);
    }
  }
  Halos.resize(it1-Halos.begin()+1);
  for(auto &&h: Halos)
    h.AverageCoordinates();
}

#include "../halo_particle_iterator.h"
static void ExchangeHalos(MpiWorker_t& world, vector <Halo_t>& InHalos, vector<Halo_t>& OutHalos, MPI_Datatype MPI_Halo_Shell_Type)
{
  typedef typename vector <Halo_t>::iterator HaloIterator_t;
  typedef HaloParticleIterator_t<HaloIterator_t> ParticleIterator_t;

  vector <IdRank_t>TargetRank(InHalos.size());
  DecideTargetProcessor(world, InHalos, TargetRank);

  //distribute halo shells
	vector <int> SendHaloCounts(world.size(),0), RecvHaloCounts(world.size()), SendHaloDisps(world.size()), RecvHaloDisps(world.size());
	sort(TargetRank.begin(), TargetRank.end(), CompareRank);
	vector <Halo_t> InHalosSorted(InHalos.size());
	vector <HBTInt> InHaloSizes(InHalos.size());
	for(HBTInt haloid=0;haloid<InHalos.size();haloid++)
	{
	  InHalosSorted[haloid]=move(InHalos[TargetRank[haloid].Id]);
	  SendHaloCounts[TargetRank[haloid].Rank]++;
	  InHaloSizes[haloid]=InHalosSorted[haloid].Particles.size();
	}
	MPI_Alltoall(SendHaloCounts.data(), 1, MPI_INT, RecvHaloCounts.data(), 1, MPI_INT, world.Communicator);
	CompileOffsets(SendHaloCounts, SendHaloDisps);
	HBTInt NumNewHalos=CompileOffsets(RecvHaloCounts, RecvHaloDisps);
	OutHalos.resize(OutHalos.size()+NumNewHalos);
	auto NewHalos=OutHalos.end()-NumNewHalos;
	MPI_Alltoallv(InHalosSorted.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_Halo_Shell_Type, &NewHalos[0], RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_Halo_Shell_Type, world.Communicator);
  //resize receivehalos
	vector <HBTInt> OutHaloSizes(NumNewHalos);
	MPI_Alltoallv(InHaloSizes.data(), SendHaloCounts.data(), SendHaloDisps.data(), MPI_HBT_INT, OutHaloSizes.data(), RecvHaloCounts.data(), RecvHaloDisps.data(), MPI_HBT_INT, world.Communicator);
	for(HBTInt i=0;i<NumNewHalos;i++)
	  NewHalos[i].Particles.resize(OutHaloSizes[i]);

	{
	//distribute halo particles
	MPI_Datatype MPI_HBT_Particle;
	Particle_t().create_MPI_type(MPI_HBT_Particle);
	//create combined iterator for each bunch of haloes
	vector <ParticleIterator_t> InParticleIterator(world.size());
	vector <ParticleIterator_t> OutParticleIterator(world.size());
	for(int rank=0;rank<world.size();rank++)
	{
	  InParticleIterator[rank].init(InHalosSorted.begin()+SendHaloDisps[rank], InHalosSorted.begin()+SendHaloDisps[rank]+SendHaloCounts[rank]);
	  OutParticleIterator[rank].init(NewHalos+RecvHaloDisps[rank], NewHalos+RecvHaloDisps[rank]+RecvHaloCounts[rank]);
	}
	vector <HBTInt> InParticleCount(world.size(),0);
	for(HBTInt i=0;i<InHalosSorted.size();i++)
	  InParticleCount[TargetRank[i].Rank]+=InHalosSorted[i].Particles.size();

	MyAllToAll<Particle_t, ParticleIterator_t, ParticleIterator_t>(world, InParticleIterator, InParticleCount, OutParticleIterator, MPI_HBT_Particle);

	MPI_Type_free(&MPI_HBT_Particle);
	}
}

void ExchangeAndMerge(MpiWorker_t& world, vector< Halo_t >& Halos)
{
  vector <Halo_t> LocalHalos;
  MPI_Datatype MPI_Halo_Shell_t;
  create_MPI_Halo_Id_type(MPI_Halo_Shell_t);
  ExchangeHalos(world, Halos, LocalHalos,  MPI_Halo_Shell_t);
  MPI_Type_free(&MPI_Halo_Shell_t);
  Halos.swap(LocalHalos);
  MergeHalos(Halos);
}
}

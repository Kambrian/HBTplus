#ifndef HALO_H_INCLUDED
#define HALO_H_INCLUDED

#include <iostream>
#include <new>
#include <algorithm>

#include "datatypes.h"
#include "snapshot_number.h"
#include "snapshot.h"
class Halo_t
{
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
	ar & Particles;
	ar & HaloId;
	ar & ComovingPosition;
	ar & PhysicalVelocity;
  }  
public:
  typedef vector <Particle_t> ParticleList_t;
  ParticleList_t Particles;
  HBTInt HaloId;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  void MoveTo(Halo_t & dest, bool MoveParticle=true)
  {
	dest.HaloId=HaloId;
	copyHBTxyz(dest.ComovingPosition, ComovingPosition);
	copyHBTxyz(dest.PhysicalVelocity, PhysicalVelocity);
	if(MoveParticle)
	  dest.Particles.swap(Particles);
  }
};

class HaloSnapshot_t: public Snapshot_t
{  
  typedef vector <Halo_t> HaloList_t;
  MPI_Datatype MPI_HBT_HaloId_t;//MPI datatype ignoring the particle list
  void BuildMPIDataType();
public:
  const ParticleSnapshot_t * ParticleSnapshot;
  HaloList_t Halos;
  HBTInt TotNumberOfParticles;
  HBTInt NumPartOfLargestHalo;
  
  HaloSnapshot_t(): Snapshot_t(), Halos(), ParticleSnapshot(nullptr), TotNumberOfParticles(0), NumPartOfLargestHalo(0)
  {
	BuildMPIDataType();
  }
  ~HaloSnapshot_t()
  {
	MPI_Type_free(&MPI_HBT_HaloId_t);
  }
  void Load(mpi::communicator & world, int snapshot_index);
  void Clear();
  void UpdateParticles(mpi::communicator & world, const ParticleSnapshot_t & snapshot);
//   void ParticleIndexToId();
  void AverageCoordinates();
  HBTInt size() const
  { 
	return Halos.size();
  }
  HBTInt GetId(HBTInt index) const
  {
	return Halos[index].HaloId;
  }
  const HBTxyz & GetComovingPosition(HBTInt index) const
  {
	return Halos[index].ComovingPosition;
  }
  const HBTxyz & GetPhysicalVelocity(HBTInt index) const
  {
	return Halos[index].PhysicalVelocity;
  }
  HBTReal GetMass(HBTInt index) const
  {
	return Halos[index].Particles.size();
  }
};

struct SizeRank_t
{
	HBTInt np;
	int rank;
};
inline bool CompareRank(const SizeRank_t &a, const SizeRank_t &b)
{
  return (a.rank<b.rank);
}
template <class Halo_T>
void DistributeHaloes(mpi::communicator &world, int root, vector <Halo_T> & InHalos, vector <Halo_T> & OutHalos, const ParticleSnapshot_t &snap, MPI_Datatype MPI_Halo_Shell_Type)
/*distribute InHalos from root to around world. 
 *the destination of each halo is the one whose particle snapshot holds the most of this halo's particles. 
 *the distributed haloes are appended to OutHalos on each node.
*/
{
  int thisrank=world.rank();
  vector <typename Halo_T::ParticleList_t> HaloBuffers;
  if(world.rank()==root)
  {
	HaloBuffers.resize(InHalos.size());
	for(HBTInt i=0;i<InHalos.size();i++)
	  HaloBuffers[i].swap(InHalos[i].Particles);
  }
  MPI_Datatype MPI_HBT_Particle;
  create_MPI_Particle_type(MPI_HBT_Particle);
 
  //broadcast haloes
  {
  HBTInt nhalo=HaloBuffers.size();
  MPI_Bcast(&nhalo, 1, MPI_HBT_INT, root, world);
  if(world.rank()!=root)
	HaloBuffers.resize(nhalo);
    
  //broadcast sizes and prepare buffer
  vector <HBTInt> HaloSizes(nhalo);
  if(world.rank()==root)
  {
	for(HBTInt i=0;i<nhalo;i++)
	  HaloSizes[i]=HaloBuffers[i].size();
  }
  MPI_Bcast(HaloSizes.data(), nhalo, MPI_HBT_INT, root, world);
  if(world.rank()!=root)
  {
	for(HBTInt i=0;i<nhalo;i++)
	  HaloBuffers[i].resize(HaloSizes[i]);
  }
  
  //broadcast particles
  vector <MPI_Aint> HaloAddress(nhalo);
  for(HBTInt i=0;i<nhalo;i++)
  {
	 MPI_Aint p;
	 MPI_Address(HaloBuffers[i].data(),&p);
	 HaloAddress[i]=p;
  }
  MPI_Datatype HaloType;
  MPI_Type_create_hindexed(nhalo, HaloSizes.data(), HaloAddress.data(), MPI_HBT_Particle, &HaloType);
  MPI_Type_commit(&HaloType);
  MPI_Bcast(MPI_BOTTOM, 1, HaloType, root, world);
  MPI_Type_free(&HaloType);
  }
  
  //fill particles and decide movement
  vector <SizeRank_t> size(HaloBuffers.size()), maxsize(HaloBuffers.size());
  for(HBTInt haloid=0;haloid<HaloBuffers.size();haloid++)
  {
	typename Halo_T::ParticleList_t & Particles=HaloBuffers[haloid];
	size[haloid].np=0;
	size[haloid].rank=thisrank;
	HBTInt &np=size[haloid].np;
	for(HBTInt i=0;i<Particles.size();i++)
	{
	  HBTInt index=snap.GetIndex(Particles[i]);
	  if(index!=SpecialConst::NullParticleId)
		Particles[np++]=snap.Particles[index];
	}
	Particles.resize(np);
  }
#ifdef HBT_INT8
#define MPI_HBTPair MPI_LONG_INT
#else
#define MPI_HBTPair MPI_2INT
#endif
	MPI_Allreduce(size.data(), maxsize.data(), size.size(), MPI_HBTPair, MPI_MAXLOC, world);
	
	vector <vector <MPI_Aint> > SendBuffers(world.size()), ReceiveBuffers(world.size());
	vector <vector<int> > SendSizes(world.size()), ReceiveSizes(world.size());
	vector <MPI_Datatype> SendTypes(world.size()), ReceiveTypes(world.size());
  
	for(HBTInt haloid=0;haloid<HaloBuffers.size();haloid++)//packing
	{
	  int rank=maxsize[haloid].rank;
	  auto & Particles=HaloBuffers[haloid];
	  MPI_Aint p;
	  MPI_Address(Particles.data(),&p);
	  SendBuffers[rank].push_back(p);
	  SendSizes[rank].push_back(Particles.size());
	}
	for(int rank=0;rank<world.size();rank++)
	{
	  if(accumulate(SendSizes[rank].begin(), SendSizes[rank].end(), 0L)>INT_MAX)
		throw runtime_error("Error: sending more than INT_MAX particles around with MPI causes overflow. try increase the number of mpi threads.\n");
	  MPI_Type_create_hindexed(SendSizes[rank].size(), SendSizes[rank].data(), SendBuffers[rank].data(), MPI_HBT_Particle, &SendTypes[rank]);
	  MPI_Type_commit(&SendTypes[rank]);
	}

	VectorAllToAll(world, SendSizes, ReceiveSizes, MPI_INT);
	
	HBTInt NumNewHalos=ReceiveSizes[0].size();
	OutHalos.resize(OutHalos.size()+NumNewHalos);
	auto NewHalos=OutHalos.end()-NumNewHalos;
	for(int rank=0;rank<world.size();rank++)
	  ReceiveBuffers[rank].resize(NumNewHalos);
	for(HBTInt haloid=0;haloid<NumNewHalos;haloid++)
	{
	  HBTInt np=0;
	  for(int rank=0; rank<world.size(); rank++)
		np+=ReceiveSizes[rank][haloid];
	  auto &Particles=NewHalos[haloid].Particles;
	  Particles.resize(np);
	  np=0;
	  for(int rank=0; rank<world.size(); rank++)
	  {
		MPI_Aint p;
		MPI_Address(Particles.data()+np, &p);
		ReceiveBuffers[rank][haloid]=p;
		np+=ReceiveSizes[rank][haloid];
	  }
	}	
	for(int rank=0;rank<world.size();rank++)
	{
	  MPI_Type_create_hindexed(NumNewHalos, ReceiveSizes[rank].data(), ReceiveBuffers[rank].data(), MPI_HBT_Particle, &ReceiveTypes[rank]);
	  MPI_Type_commit(&ReceiveTypes[rank]);
	}
	
	vector <int> Counts(world.size(),1), Disps(world.size(),0);
	MPI_Alltoallw(MPI_BOTTOM, Counts.data(), Disps.data(), SendTypes.data(), MPI_BOTTOM, Counts.data(), Disps.data(), ReceiveTypes.data(), world);
	for(int rank=0;rank<world.size();rank++)
	{
	  MPI_Type_free(&SendTypes[rank]);
	  MPI_Type_free(&ReceiveTypes[rank]);
	}
	MPI_Type_free(&MPI_HBT_Particle);
	
	//copy other properties
	vector <Halo_T> TmpHalos;
	if(world.rank()==root)
	{//reuse maxsize for sorting
	  for(HBTInt haloid=0;haloid<InHalos.size();haloid++)
		maxsize[haloid].np=haloid;
	  stable_sort(maxsize.begin(), maxsize.end(), CompareRank);
	  TmpHalos.resize(InHalos.size());
	  for(HBTInt haloid=0;haloid<InHalos.size();haloid++)
		InHalos[maxsize[haloid].np].MoveTo(TmpHalos[haloid], false);
	  for(int rank=0;rank<world.size();rank++)
		Counts[rank]=SendSizes[rank].size();
	  CompileOffsets(Counts, Disps);
	}
	MPI_Scatterv(TmpHalos.data(), Counts.data(), Disps.data(), MPI_Halo_Shell_Type, &NewHalos[0], NumNewHalos, MPI_Halo_Shell_Type, root, world);
}
#endif
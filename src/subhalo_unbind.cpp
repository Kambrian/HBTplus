#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"
#include "gravity_tree.h"

struct ParticleEnergy_t
{
  HBTInt pid;
  float E;
};
inline bool CompEnergy(const ParticleEnergy_t & a, const ParticleEnergy_t & b)
{
  return (a.E<b.E);
};
static HBTInt PartitionBindingEnergy(vector <ParticleEnergy_t> &Elist, const size_t len)
/*sort Elist to move unbound particles to the end*/
{//similar to the C++ partition() func
  if(len==0) return 0;
  if(len==1) return Elist[0].E<0;
  
  ParticleEnergy_t Etmp=Elist[0]; 
  auto iterforward=Elist.begin(), iterbackward=Elist.begin()+len;
  while(true)
  {
	//iterforward is a void now, can be filled
	while(true)
	{
	  iterbackward--;
	  if(iterbackward==iterforward)
	  {
		*iterforward=Etmp;
		if(Etmp.E<0) iterbackward++;
		return iterbackward-Elist.begin();
	  }
	  if(iterbackward->E<0) break;
	}
	*iterforward=*iterbackward;
	//iterbackward is a void now, can be filled
	while(true)
	{
	  iterforward++;
	  if(iterforward==iterbackward)
	  {
		*iterbackward=Etmp;
		if(Etmp.E<0) iterbackward++;
		return iterbackward-Elist.begin();
	  }
	  if(iterforward->E>0) break;
	}
	*iterbackward=*iterforward;
  }
}
static void PopMostBoundParticle(ParticleEnergy_t * Edata, const HBTInt Nbound)
{
  HBTInt imin=0;
  for(HBTInt i=1;i<Nbound;i++)
  {
	if(Edata[i].E<Edata[imin].E) imin=i;
  }
  if(imin!=0) swap(Edata[imin], Edata[0]);
}
class EnergySnapshot_t: public Snapshot_t
{
  HBTInt GetParticle(HBTInt i) const
  {
	return Elist[i].pid;
  }
public:
  ParticleEnergy_t * Elist;
  typedef vector <Particle_t> ParticleList_t;
  HBTInt N;
  const ParticleList_t & Particles;
  EnergySnapshot_t(ParticleEnergy_t *e, HBTInt n, const ParticleList_t &particles, const Snapshot_t & epoch): Elist(e), N(n), Particles(particles)
  {
	SetEpoch(epoch);
  };
  HBTInt size() const
  {
	return N;
  }
  HBTInt GetId(HBTInt i) const
  {
	return Particles[GetParticle(i)].Id;
  }
  HBTReal GetMass(HBTInt i) const
  {
	return Particles[GetParticle(i)].Mass;
  }
  const HBTxyz & GetPhysicalVelocity(HBTInt i) const
  {
	return Particles[GetParticle(i)].PhysicalVelocity;
  }
  const HBTxyz & GetComovingPosition(HBTInt i) const
  {
	return Particles[GetParticle(i)].ComovingPosition;
  }
};
void Subhalo_t::Unbind(const Snapshot_t &epoch)
{//the reference frame (pos and vel) should already be initialized before unbinding.
  if(1==Particles.size()) return;
  
  const HBTInt OldMostBoundParticle=0;
  OctTree_t tree;
  tree.Reserve(Particles.size());
  Nbound=Particles.size(); //start from full set
  HBTInt Nlast; 
  
  vector <ParticleEnergy_t> Elist(Nbound);
	for(HBTInt i=0;i<Nbound;i++)
	  Elist[i].pid=i;
  EnergySnapshot_t ESnap(Elist.data(), Elist.size(), Particles, epoch);
	while(true)
	{
		Nlast=Nbound;
		tree.Build(ESnap, Nlast);
	  #pragma omp parallel for if(Nlast>100)
	  for(HBTInt i=0;i<Nlast;i++)
	  {
		HBTInt pid=Elist[i].pid;
		Elist[i].E=tree.BindingEnergy(Particles[pid].ComovingPosition, Particles[pid].PhysicalVelocity, ComovingPosition, PhysicalVelocity, Particles[pid].Mass);
	  }
		Nbound=PartitionBindingEnergy(Elist, Nlast);//TODO: parallelize this.
		if(Nbound<HBTConfig.MinNumPartOfSub)
		{
		  Nbound=0;
		  Elist[0].pid=OldMostBoundParticle;
		}
		else
		{
		  sort(Elist.begin()+Nbound, Elist.begin()+Nlast, CompEnergy); //only sort the unbound part
		  PopMostBoundParticle(Elist.data(), Nbound);
		}
		copyHBTxyz(ComovingPosition, Particles[Elist[0].pid].ComovingPosition);
		copyHBTxyz(PhysicalVelocity, Particles[Elist[0].pid].PhysicalVelocity);
		if(0==Nbound||Nbound>Nlast*HBTConfig.BoundMassPrecision)  break;
	}
	if(Nbound)
	{
	  sort(Elist.begin(), Elist.begin()+Nbound, CompEnergy); //sort the self-bound part
	  Nlast=Nbound*HBTConfig.SourceSubRelaxFactor;
	  if(Nlast>Particles.size()) Nlast=Particles.size();
	  //todo: optimize this with in-place permutation, to avoid mem alloc and copying.
	  ParticleList_t p(Nlast);
	  for(HBTInt i=0;i<Nlast;i++) 
		p[i]=Particles[Elist[i].pid];
	  Particles.swap(p);
	}
	else
	{
	  Nbound=1;
	  Nlast=1;//what if this is a central?? any fixes?
	  Particles.resize(Nlast);//keep old particle order, just shrink.
	}
  //todo: output angular momentum and total energy as well, for calculation of spin.
}
void SubhaloSnapshot_t::RefineParticles()
{//it's more expensive to build an exclusive list. so do inclusive here. 
  //TODO: ensure the inclusive unbinding is stable (contaminating particles from big subhaloes may hurdle the unbinding
#ifdef _OPENMP
 if(ParallelizeHaloes) cout<<"Unbinding with HaloPara...\n";
 else cout<<"Unbinding with ParticlePara...\n";
#else
 cout<<"Unbinding..."<<endl;
#endif  
#pragma omp parallel for if(ParallelizeHaloes)
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	try
	{
	  Subhalos[subid].Unbind(*this);
	}
	catch(OctTreeExceeded_t &tree_exception)
	{
	  cerr<<"Error: "<<tree_exception.what()<<" in subhalo "<<subid<<endl;
	  exit(1);
	}
  }
}
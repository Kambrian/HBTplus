//TODO: unify the reference frame for specificProperties...
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
  void AverageVelocity(HBTxyz& CoV, HBTInt NumPart)
  /*mass weighted average velocity*/
  {
	HBTInt i,j;
	double sv[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoV, GetPhysicalVelocity(0));
	  return;
	}
	
	sv[0]=sv[1]=sv[2]=0.;
	msum=0.;
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=GetMass(i);
	  msum+=m;
	  for(j=0;j<3;j++)
		sv[j]+=GetPhysicalVelocity(i)[j]*m;
	}
	
	for(j=0;j<3;j++)
	  CoV[j]=sv[j]/msum;
  }
  void AveragePosition(HBTxyz& CoM, HBTInt NumPart)
  /*mass weighted average position*/
  {
	HBTInt i,j;
	double sx[3],origin[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoM, GetComovingPosition(0));
	  return;
	}
	
	sx[0]=sx[1]=sx[2]=0.;
	msum=0.;
	if(HBTConfig.PeriodicBoundaryOn)
	  for(j=0;j<3;j++)
		origin[j]=GetComovingPosition(0)[j];
	  
	  for(i=0;i<NumPart;i++)
	  {
		HBTReal m=GetMass(i);
		msum+=m;
		for(j=0;j<3;j++)
		  if(HBTConfig.PeriodicBoundaryOn)
			sx[j]+=NEAREST(GetComovingPosition(i)[j]-origin[j])*m;
		  else
			sx[j]+=GetComovingPosition(i)[j]*m;
	  }
	  
	  for(j=0;j<3;j++)
	  {
		sx[j]/=msum;
		if(HBTConfig.PeriodicBoundaryOn) sx[j]+=origin[j];
		CoM[j]=sx[j];
	  }
  }
  void AverageKinematics(float &SpecificPotentialEnergy, float &SpecificKineticEnergy, float SpecificAngularMomentum[3], HBTInt NumPart, const HBTxyz & refPos, const HBTxyz &refVel)
  /*obtain specific potential, kinetic energy, and angular momentum for the first NumPart particles
   * currently only support fixed particle mass
   * 
   * TODO: extend to variable mass?
   * Note there is a slight inconsistency in the energy since they were calculated from the previous unbinding loop, but the refVel has been updated.
   */
  {
	if(NumPart<=1)
	{
	  SpecificPotentialEnergy=0.;
	  SpecificKineticEnergy=0.;
	  SpecificAngularMomentum[0]=SpecificAngularMomentum[1]=SpecificAngularMomentum[2]=0.;
	  return;
	}
	double E=0., K=0., AM[3]={0.};
	for(HBTInt i=0;i<NumPart;i++)
	{
	  E+=Elist[i].E;
	  const HBTxyz & x=GetComovingPosition(i);
	  const HBTxyz & v=GetPhysicalVelocity(i);
	  double dx[3], dv[3];
	  for(int j=0;j<3;j++)
	  {
		dx[j]=x[j]-refPos[j];
		if(HBTConfig.PeriodicBoundaryOn) dx[j]=NEAREST(dx[j]);
		dx[j]*=ScaleFactor; //physical
		dv[j]=v[j]-refVel[j]+Hz*dx[j];
		K+=dv[j]*dv[j];
	  }
	  AM[0]+=dx[1]*dv[2]-dx[2]*dv[1];
	  AM[1]+=dx[0]*dv[2]-dx[2]*dv[0];
	  AM[2]+=dx[0]*dv[1]-dx[1]*dv[0];
	}
	E/=NumPart;
	K*=0.5/NumPart;
	SpecificPotentialEnergy=E-K;
	SpecificKineticEnergy=K;
	for(int j=0;j<3;j++) 
	  SpecificAngularMomentum[j]=AM[j]/NumPart;	//physical
  }
};
void Subhalo_t::Unbind(const Snapshot_t &epoch)
{//the reference frame (pos and vel) should already be initialized before unbinding.
  if(1==Particles.size()) return;
  
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
		Elist[i].E=tree.BindingEnergy(Particles[pid].ComovingPosition, Particles[pid].PhysicalVelocity, ComovingMostBoundPosition, PhysicalAverageVelocity, Particles[pid].Mass);
	  }
		Nbound=PartitionBindingEnergy(Elist, Nlast);//TODO: parallelize this.
		if(Nbound<HBTConfig.MinNumPartOfSub)//disruption
		{
		  Nbound=1;
		  Nlast=1;
		  Particles.resize(1);//old particle list retained, only need to shrink
		  SnapshotIndexOfDeath=epoch.GetSnapshotIndex();
		  copyHBTxyz(ComovingMostBoundPosition, Particles[0].ComovingPosition);
		  copyHBTxyz(PhysicalMostBoundVelocity, Particles[0].PhysicalVelocity);
		  copyHBTxyz(ComovingAveragePosition, ComovingMostBoundPosition);
		  copyHBTxyz(PhysicalAverageVelocity, PhysicalMostBoundVelocity);
		  break;
		}
		else
		{
		  sort(Elist.begin()+Nbound, Elist.begin()+Nlast, CompEnergy); //only sort the unbound part
		  PopMostBoundParticle(Elist.data(), Nbound);
		  ESnap.AverageVelocity(PhysicalAverageVelocity, Nbound);
		  copyHBTxyz(ComovingMostBoundPosition, Particles[Elist[0].pid].ComovingPosition);
		  if(Nbound>Nlast*HBTConfig.BoundMassPrecision)//converge
		  {
			//update properties
			ESnap.AveragePosition(ComovingAveragePosition, Nbound);
			copyHBTxyz(PhysicalMostBoundVelocity, Particles[Elist[0].pid].PhysicalVelocity);
			//update particle list
			sort(Elist.begin(), Elist.begin()+Nbound, CompEnergy); //sort the self-bound part
			Nlast=Nbound*HBTConfig.SourceSubRelaxFactor;
			if(Nlast>Particles.size()) Nlast=Particles.size();
			//todo: optimize this with in-place permutation, to avoid mem alloc and copying.
			ParticleList_t p(Nlast);
			for(HBTInt i=0;i<Nlast;i++)
			{
			  p[i]=Particles[Elist[i].pid];
			  Elist[i].pid=i;//update particle index in Elist as well.
			}
			Particles.swap(p);
			break;
		  }
		}
	}
	ESnap.AverageKinematics(SpecificSelfPotentialEnergy, SpecificSelfKineticEnergy, SpecificAngularMomentum, Nbound, ComovingMostBoundPosition, PhysicalAverageVelocity);
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
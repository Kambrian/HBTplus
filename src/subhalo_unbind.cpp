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
	double svx,svy,svz,msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoV, GetPhysicalVelocity(0));
	  return;
	}
	
	svx=svy=svz=0.;
	msum=0.;
	#pragma omp paralle for reduction(+:msum, svx, svy, svz) if(NumPart>100)
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=GetMass(i);
	  const HBTxyz &v=GetPhysicalVelocity(i);
	  msum+=m;
	  svx+=v[0]*m;
	  svy+=v[1]*m;
	  svz+=v[2]*m;
	}
	
	CoV[0]=svx/msum;
	CoV[1]=svy/msum;
	CoV[2]=svz/msum;
  }
  void AveragePosition(HBTxyz& CoM, HBTInt NumPart)
  /*mass weighted average position*/
  {
	HBTInt i,j;
	double sx,sy,sz,origin[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoM, GetComovingPosition(0));
	  return;
	}
	
	if(HBTConfig.PeriodicBoundaryOn)
	  for(j=0;j<3;j++)
		origin[j]=GetComovingPosition(0)[j];
	
	sx=sy=sz=0.;
	msum=0.;
	#pragma omp paralle for reduction(+:msum, sx, sy, sz) if(NumPart>100)
	  for(i=0;i<NumPart;i++)
	  {
		HBTReal m=GetMass(i);
		const HBTxyz &x=GetComovingPosition(i);
		msum+=m;
		if(HBTConfig.PeriodicBoundaryOn)
		{
		  sx+=NEAREST(x[0]-origin[0])*m;
		  sy+=NEAREST(x[1]-origin[1])*m;
		  sz+=NEAREST(x[2]-origin[2])*m;
		}
		else
		{
		  sx+=x[0]*m;
		  sy+=x[1]*m;
		  sz+=x[2]*m;
		}
	  }
	  sx/=msum; sy/=msum; sz/=msum;
	  if(HBTConfig.PeriodicBoundaryOn)
	  {
		sx+=origin[0];
		sy+=origin[1];
		sz+=origin[2];
	  }
	  CoM[0]=sx;
	  CoM[1]=sy;
	  CoM[2]=sz;
  }
  void AverageKinematics(float &SpecificPotentialEnergy, float &SpecificKineticEnergy, float SpecificAngularMomentum[3], HBTInt NumPart, const HBTxyz & refPos, const HBTxyz &refVel)
  /*obtain specific potential, kinetic energy, and angular momentum for the first NumPart particles
   * all quantities are physical
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
	double E=0., K=0., AMx=0., AMy=0., AMz=0.;
	#pragma omp parallel for reduction(+:E, K, AMx, AMy, AMz) if(NumPart>100)
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
	  AMx+=dx[1]*dv[2]-dx[2]*dv[1];
	  AMy+=dx[0]*dv[2]-dx[2]*dv[0];
	  AMz+=dx[0]*dv[1]-dx[1]*dv[0];
	}
	E/=NumPart;
	K*=0.5/NumPart;
	SpecificPotentialEnergy=E-K;
	SpecificKineticEnergy=K;
	SpecificAngularMomentum[0]=AMx/NumPart;
	SpecificAngularMomentum[1]=AMy/NumPart;
	SpecificAngularMomentum[2]=AMz/NumPart;
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
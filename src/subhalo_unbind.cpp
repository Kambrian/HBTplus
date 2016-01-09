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
public:
  ParticleEnergy_t * Elist;
  HBTInt N;
  const Snapshot_t & Snapshot;
  EnergySnapshot_t(ParticleEnergy_t *e, HBTInt n, const Snapshot_t & fullsnapshot): Elist(e), N(n), Snapshot(fullsnapshot)
  {
	SetEpoch(fullsnapshot);
  };
  HBTInt size() const
  {
	return N;
  }
  HBTInt GetMemberId(const HBTInt i) const
  {
	return Elist[i].pid;
  }
  HBTReal GetMass(const HBTInt i) const
  {
	return Snapshot.GetMass(GetMemberId(i));
  }
  const HBTxyz & GetPhysicalVelocity(const HBTInt i) const
  {
	return Snapshot.GetPhysicalVelocity(GetMemberId(i));
  }
  const HBTxyz & GetComovingPosition(const HBTInt i) const
  {
	return Snapshot.GetComovingPosition(GetMemberId(i));
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
  void AveragePosition(HBTxyz& CoM, HBTInt NumPart) const
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
  void AverageKinematics(HBTReal &SpecificPotentialEnergy, HBTReal &SpecificKineticEnergy, HBTxyz &SpecificAngularMomentum, HBTInt NumPart, const HBTxyz & refPos, const HBTxyz &refVel)
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
void Subhalo_t::Unbind(const ParticleSnapshot_t &snapshot)
{//the reference frame should already be initialized before unbinding.
  if(1==Particles.size()) return;
  
  HBTInt OldMostBoundParticle=Particles[0];
  OctTree_t tree;
  tree.Reserve(Particles.size());
  Nbound=Particles.size(); //start from full set
  HBTInt Nlast; 
  
  vector <ParticleEnergy_t> Elist(Nbound);
	for(HBTInt i=0;i<Nbound;i++)
	  Elist[i].pid=Particles[i];
  EnergySnapshot_t ESnap(Elist.data(), Elist.size(), snapshot);
	while(true)
	{
		Nlast=Nbound;
		tree.Build(ESnap, Nlast);
	  #pragma omp parallel for if(Nlast>100)
	  for(HBTInt i=0;i<Nlast;i++)
	  {
		HBTInt pid=Elist[i].pid;
		Elist[i].E=tree.BindingEnergy(snapshot.GetComovingPosition(pid), snapshot.GetPhysicalVelocity(pid), ComovingMostBoundPosition, PhysicalAverageVelocity, snapshot.GetParticleMass(pid));
	  }
		Nbound=PartitionBindingEnergy(Elist, Nlast);//TODO: parallelize this.
		if(Nbound<HBTConfig.MinNumPartOfSub)//disruption
		{
		  Nbound=1;
		  SnapshotIndexOfDeath=snapshot.GetSnapshotIndex();
		  Particles[0]=OldMostBoundParticle;
		  Particles.resize(1);//what if this is a central?? downgrade to sat later.
		  copyHBTxyz(ComovingMostBoundPosition, snapshot.GetComovingPosition(OldMostBoundParticle));
		  copyHBTxyz(PhysicalMostBoundVelocity, snapshot.GetPhysicalVelocity(OldMostBoundParticle));
		  copyHBTxyz(ComovingAveragePosition, ComovingMostBoundPosition);
		  copyHBTxyz(PhysicalAverageVelocity, PhysicalMostBoundVelocity);
		  break;
		}
		else
		{
		  sort(Elist.begin()+Nbound, Elist.begin()+Nlast, CompEnergy); //only sort the unbound part
		  PopMostBoundParticle(Elist.data(), Nbound);
		  ESnap.AverageVelocity(PhysicalAverageVelocity, Nbound);
		  copyHBTxyz(ComovingMostBoundPosition, snapshot.GetComovingPosition(Elist[0].pid));
		  if(Nbound>Nlast*HBTConfig.BoundMassPrecision)//converge
		  {
			sort(Elist.begin(), Elist.begin()+Nbound, CompEnergy); //sort the self-bound part
			//update particle list
			Nlast=Nbound*HBTConfig.SourceSubRelaxFactor;
			if(Nlast>Particles.size()) Nlast=Particles.size();
			for(HBTInt i=0;i<Nlast;i++) Particles[i]=Elist[i].pid;
			Particles.resize(Nlast);
			//update properties
			ESnap.AveragePosition(ComovingAveragePosition, Nbound);
			copyHBTxyz(PhysicalMostBoundVelocity, snapshot.GetPhysicalVelocity(Elist[0].pid));
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
	  Subhalos[subid].Unbind(*SnapshotPointer);
	}
	catch(OctTreeExceeded_t &tree_exception)
	{
	  cerr<<"Error: "<<tree_exception.what()<<" in subhalo "<<subid<<endl;
	  exit(1);
	}
  }
}
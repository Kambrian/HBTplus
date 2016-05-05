//TODO: unify the reference frame for specificProperties...
#include <iostream>
#include <new>
#include <algorithm>
#include <boost/concept_check.hpp>
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
  HBTReal MassFactor;
  EnergySnapshot_t(ParticleEnergy_t *e, HBTInt n, const Snapshot_t & fullsnapshot): Elist(e), N(n), Snapshot(fullsnapshot), MassFactor(1.)
  {
	SetEpoch(fullsnapshot);
  };
  void SetMassUnit(HBTReal mass_unit)
  {
	MassFactor=mass_unit;
  }
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
	return Snapshot.GetMass(GetMemberId(i))*MassFactor;
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
  void AveragePosition(HBTxyz& CoM, HBTInt NumPart) const
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
inline void RefineBindingEnergyOrder(EnergySnapshot_t &ESnap, HBTInt Size, OctTree_t &tree, HBTxyz &RefPos, HBTxyz &RefVel)
{//reorder the first Size particles according to their self-binding energy
  auto &Elist=ESnap.Elist;
  auto &snapshot=ESnap.Snapshot;
  tree.Build(ESnap, Size);
  vector <ParticleEnergy_t> Einner(Size);
  #pragma omp parallel if(Size>100)
  {
	#pragma omp for
	for(HBTInt i=0;i<Size;i++)
	{
	  HBTInt pid=Elist[i].pid;
	  Einner[i].pid=i;
	  Einner[i].E=tree.BindingEnergy(snapshot.GetComovingPosition(pid), snapshot.GetPhysicalVelocity(pid), RefPos, RefVel, snapshot.GetMass(pid));
	}
	#pragma omp single
	sort(Einner.begin(), Einner.end(), CompEnergy);
	#pragma omp for
	for(HBTInt i=0;i<Size;i++)
	{
	  Einner[i]=Elist[Einner[i].pid];
	}
	#pragma omp for
	for(HBTInt i=0;i<Size;i++)
	{
	  Elist[i]=Einner[i];
	}
  }
}

void Subhalo_t::Unbind(const ParticleSnapshot_t &snapshot)
{//the reference frame should already be initialized before unbinding.
  HBTInt MaxSampleSize=HBTConfig.MaxSampleSizeOfPotentialEstimate;
  bool RefineMostboundParticle=(MaxSampleSize>0&&HBTConfig.RefineMostboundParticle);
  HBTReal BoundMassPrecision=HBTConfig.BoundMassPrecision;
  
  if(1==Particles.size()) return;
  HBTxyz OldRefPos, OldRefVel;
  auto &RefPos=ComovingAveragePosition;
  auto &RefVel=PhysicalAverageVelocity;
  
  OctTree_t tree;
  tree.Reserve(Particles.size());
  Nbound=Particles.size(); //start from full set
  random_shuffle(Particles.begin(), Particles.end()); //shuffle for easy resampling later.
  HBTInt Nlast; 
  
  vector <ParticleEnergy_t> Elist(Nbound);
	for(HBTInt i=0;i<Nbound;i++)
	  Elist[i].pid=Particles[i];
  EnergySnapshot_t ESnap(Elist.data(), Elist.size(), snapshot);
  bool CorrectionLoop=false;
	while(true)
	{
	  if(CorrectionLoop)
	  {//correct the potential due to removed particles
		#define VecNorm(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
		#define VecDot(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
		HBTxyz RefVelDiff;
		snapshot.RelativeVelocity(OldRefPos, OldRefVel, RefPos, RefVel, RefVelDiff);
		HBTReal dK=0.5*VecNorm(RefVelDiff);
		EnergySnapshot_t ESnapCorrection(&Elist[Nbound], Nlast-Nbound, snapshot); //point to freshly removed particles
		tree.Build(ESnapCorrection); 
		#pragma omp parallel for if(Nlast>100)
		for(HBTInt i=0;i<Nbound;i++)
		{
		  HBTInt pid=Elist[i].pid;
		  auto &x=snapshot.GetComovingPosition(pid);
		  auto &v=snapshot.GetPhysicalVelocity(pid);
		  HBTxyz OldVel;
		  snapshot.RelativeVelocity(x,v,OldRefPos, OldRefVel, OldVel);
		  Elist[i].E+=VecDot(OldVel, RefVelDiff)+dK-tree.EvaluatePotential(x, 0);;
		}
		#undef VecNorm
		#undef VecDot
		Nlast=Nbound;
	  }
	  else
	  {
		Nlast=Nbound;
		HBTInt np_tree=Nlast;
		if(MaxSampleSize>0&&Nlast>MaxSampleSize)//downsample
		{
		  np_tree=MaxSampleSize;
		  ESnap.SetMassUnit((HBTReal)Nlast/MaxSampleSize);
		}
		tree.Build(ESnap, np_tree);
		#pragma omp parallel for if(Nlast>100)
		for(HBTInt i=0;i<Nlast;i++)
		{
		  HBTInt pid=Elist[i].pid;
		  HBTReal mass;
		  if(i<MaxSampleSize)
			mass=ESnap.GetMass(i); //to correct for self-gravity
		  else
			mass=0.;//not sampled in tree, no self gravity to correct
		  Elist[i].E=tree.BindingEnergy(snapshot.GetComovingPosition(pid), snapshot.GetPhysicalVelocity(pid), RefPos, RefVel, mass);
		}
		ESnap.SetMassUnit(1.);//reset, no necessary
	  }
		Nbound=PartitionBindingEnergy(Elist, Nlast);//TODO: parallelize this.
		if(Nbound<HBTConfig.MinNumPartOfSub)//disruption
		{
		  Nbound=1;
		  Nlast=1;
		  SnapshotIndexOfDeath=snapshot.GetSnapshotIndex();
		  //old particle list retained. old mostbound coordinates also retained.
		  copyHBTxyz(ComovingAveragePosition, ComovingMostBoundPosition);
		  copyHBTxyz(PhysicalAverageVelocity, PhysicalMostBoundVelocity);
		  break;
		}
		else
		{
		  sort(Elist.begin()+Nbound, Elist.begin()+Nlast, CompEnergy); //only sort the unbound part
		  HBTInt Ndiff=Nlast-Nbound;
		  if(Ndiff<Nbound)
		  {
			if(MaxSampleSize<=0||Ndiff<MaxSampleSize)
			{
			  CorrectionLoop=true;
			  copyHBTxyz(OldRefPos, RefPos);
			  copyHBTxyz(OldRefVel, RefVel);
			}
		  }
		  ESnap.AverageVelocity(PhysicalAverageVelocity, Nbound);
		  ESnap.AveragePosition(ComovingAveragePosition, Nbound);
		  if(Nbound>Nlast*HBTConfig.BoundMassPrecision)//converge
		  {
			sort(Elist.begin(), Elist.begin()+Nbound, CompEnergy); //sort the self-bound part
			if(RefineMostboundParticle&&Nbound>MaxSampleSize)//refine most-bound particle, not necessary usually..
			  RefineBindingEnergyOrder(ESnap, MaxSampleSize, tree, RefPos, RefVel);
			//update particle list
			Nlast=Nbound*HBTConfig.SourceSubRelaxFactor;
			if(Nlast>Particles.size()) Nlast=Particles.size();
			for(HBTInt i=0;i<Nlast;i++) Particles[i]=Elist[i].pid;
			//update mostbound coordinate
			copyHBTxyz(ComovingMostBoundPosition, snapshot.GetComovingPosition(Elist[0].pid));
			copyHBTxyz(PhysicalMostBoundVelocity, snapshot.GetPhysicalVelocity(Elist[0].pid));
			break;
		  }
		}
	}
	Particles.resize(Nlast);
	ESnap.AverageKinematics(SpecificSelfPotentialEnergy, SpecificSelfKineticEnergy, SpecificAngularMomentum, Nbound, RefPos, RefVel);//only use CoM frame when unbinding and calculating Kinematics
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
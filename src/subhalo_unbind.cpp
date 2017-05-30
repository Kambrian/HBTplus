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
  HBTReal MassFactor;
  EnergySnapshot_t(ParticleEnergy_t *e, HBTInt n, const Snapshot_t & fullsnapshot): Elist(e), N(n), Snapshot(fullsnapshot), MassFactor(1.)
  {
	Cosmology=fullsnapshot.Cosmology;
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
  double AverageVelocity(HBTxyz& CoV, HBTInt NumPart)
  /*mass weighted average velocity*/
  {
	HBTInt i,j;
	double svx,svy,svz,msum;
	
	if(0==NumPart) return 0.;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoV, GetPhysicalVelocity(0));
	  return GetMass(0);
	}
	
	svx=svy=svz=0.;
	msum=0.;
	#pragma omp parallel for reduction(+:msum, svx, svy, svz) if(NumPart>100)
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
	return msum;
  }
  double AveragePosition(HBTxyz& CoM, HBTInt NumPart) const
  /*mass weighted average position*/
  {
	HBTInt i,j;
	double sx,sy,sz,origin[3],msum;
	
	if(0==NumPart) return 0.;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoM, GetComovingPosition(0));
	  return GetMass(0);
	}
	
	if(HBTConfig.PeriodicBoundaryOn)
	  for(j=0;j<3;j++)
		origin[j]=GetComovingPosition(0)[j];
	
	sx=sy=sz=0.;
	msum=0.;
	#pragma omp parallel for reduction(+:msum, sx, sy, sz) if(NumPart>100)
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
	  return msum;
  }
  void AverageKinematics(float &SpecificPotentialEnergy, float &SpecificKineticEnergy, float SpecificAngularMomentum[3], HBTInt NumPart, const HBTxyz & refPos, const HBTxyz &refVel)
  /*obtain specific potential, kinetic energy, and angular momentum for the first NumPart particles
   * all quantities are physical
   * 
   * Note there is a slight inconsistency in the energy since they were calculated from the previous unbinding loop, but the refVel has been updated.
   */
  {
	if(NumPart<=1)
	{
// 	  if(NumPart==1) Mbound=GetMass(0);
// 	  else Mbound=0.;
	  SpecificPotentialEnergy=0.;
	  SpecificKineticEnergy=0.;
	  SpecificAngularMomentum[0]=SpecificAngularMomentum[1]=SpecificAngularMomentum[2]=0.;
	  return;
	}
	double E=0., K=0., AMx=0., AMy=0., AMz=0., M=0.;
	#pragma omp parallel for reduction(+:E, K, AMx, AMy, AMz, M) if(NumPart>100)
	for(HBTInt i=0;i<NumPart;i++)
	{
	  HBTReal m=GetMass(i);
	  E+=Elist[i].E*m;
	  const HBTxyz & x=GetComovingPosition(i);
	  const HBTxyz & v=GetPhysicalVelocity(i);
	  double dx[3], dv[3];
	  for(int j=0;j<3;j++)
	  {
		dx[j]=x[j]-refPos[j];
		if(HBTConfig.PeriodicBoundaryOn) dx[j]=NEAREST(dx[j]);
		dx[j]*=Cosmology.ScaleFactor; //physical
		dv[j]=v[j]-refVel[j]+Cosmology.Hz*dx[j];
		K+=dv[j]*dv[j]*m;
	  }
	  AMx+=(dx[1]*dv[2]-dx[2]*dv[1])*m;
	  AMy+=(dx[0]*dv[2]-dx[2]*dv[0])*m;
	  AMz+=(dx[0]*dv[1]-dx[1]*dv[0])*m;
	  M+=m;
	}
	E/=M;
	K*=0.5/M;
	SpecificPotentialEnergy=E-K;
	SpecificKineticEnergy=K;
	SpecificAngularMomentum[0]=AMx/M;
	SpecificAngularMomentum[1]=AMy/M;
	SpecificAngularMomentum[2]=AMz/M;
// 	Mbound=M;
  }
};
inline void RefineBindingEnergyOrder(EnergySnapshot_t &ESnap, HBTInt Size, GravityTree_t &tree, HBTxyz &RefPos, HBTxyz &RefVel)
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
	  #ifdef UNBIND_WITH_THERMAL_ENERGY
	  Elist[i].E+=snapshot.GetInternalEnergy(pid);
	  #endif
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
/* remove unbound particles.
 * this will update Nbound, Mbound, NboundType, MboundType, the particle list, and kinematic properties (specific potential, kinetic energy, angular momentum, the most bound and average positions and velocities)
 * the returned particle list will be ordered by binding energy.
 * 
 * the reference frame should already be initialized before unbinding.
 */
{
  HBTInt MaxSampleSize=HBTConfig.MaxSampleSizeOfPotentialEstimate;
  bool RefineMostboundParticle=(MaxSampleSize>0&&HBTConfig.RefineMostboundParticle);
  HBTReal BoundMassPrecision=HBTConfig.BoundMassPrecision;
  
  if(Particles.size()==0)
  {
    Nbound=0;
    Mbound=0.;
#ifdef SAVE_BINDING_ENERGY
	Energies.clear();
#endif
    return;
  }
  if(Particles.size()==1) 
  {
	Nbound=1;
	Mbound=snapshot.GetMass(Particles[0]);
#ifdef SAVE_BINDING_ENERGY
	Energies.resize(1);
	Energies[0]=0.;
#endif
	return;
  }

  HBTxyz OldRefPos, OldRefVel;
  auto &RefPos=ComovingAveragePosition;
  auto &RefVel=PhysicalAverageVelocity;
  
  HBTInt OldMostboundParticle=Particles[0];//backup 
  GravityTree_t tree;
  tree.Reserve(Particles.size());
  Nbound=Particles.size(); //start from full set
  if(MaxSampleSize>0&&Nbound>MaxSampleSize)	random_shuffle(Particles.begin(), Particles.end()); //shuffle for easy resampling later.
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
		  Elist[i].E+=VecDot(OldVel, RefVelDiff)+dK-tree.EvaluatePotential(x, 0);
		}
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
		  if(i<np_tree)
			mass=ESnap.GetMass(i); //to correct for self-gravity
		  else
			mass=0.;//not sampled in tree, no self gravity to correct
		  Elist[i].E=tree.BindingEnergy(snapshot.GetComovingPosition(pid), snapshot.GetPhysicalVelocity(pid), RefPos, RefVel, mass);
#ifdef UNBIND_WITH_THERMAL_ENERGY
		  Elist[i].E+=snapshot.GetInternalEnergy(pid);
#endif
		}
		ESnap.SetMassUnit(1.);//reset
	  }
		Nbound=PartitionBindingEnergy(Elist, Nlast);//TODO: parallelize this.
		if(Nbound<HBTConfig.MinNumPartOfSub)//disruption
		{
		  Nbound=1;
		  Nlast=1;
		  if(IsAlive())
		    SnapshotIndexOfDeath=snapshot.GetSnapshotIndex();
		  for(auto &&p: Particles)
		  {
			if(p==OldMostboundParticle)
			{
			  swap(p, Particles[0]);//restore old most-bound to beginning
			  break;
			}
		  }
		  //old mostbound coordinates retained
		  copyHBTxyz(ComovingAveragePosition, ComovingMostBoundPosition);
		  copyHBTxyz(PhysicalAverageVelocity, PhysicalMostBoundVelocity);
		  Mbound=snapshot.GetMass(Particles[0]);
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
		  Mbound=ESnap.AverageVelocity(PhysicalAverageVelocity, Nbound);
		  ESnap.AveragePosition(ComovingAveragePosition, Nbound);
		  if(Nbound>=Nlast*BoundMassPrecision)//converge
		  {
			if(!IsAlive())
			  SnapshotIndexOfDeath=SpecialConst::NullSnapshotId;//clear death snapshot
			if(IsTrapped())
			  SinkTrackId=SpecialConst::NullTrackId;//clear sinktrack as well
			sort(Elist.begin(), Elist.begin()+Nbound, CompEnergy); //sort the self-bound part
			if(RefineMostboundParticle&&Nbound>MaxSampleSize)//refine most-bound particle, not necessary usually..
			  RefineBindingEnergyOrder(ESnap, MaxSampleSize, tree, RefPos, RefVel);
			//update particle list
			for(HBTInt i=0;i<Particles.size();i++) Particles[i]=Elist[i].pid;
			//update mostbound coordinate
			copyHBTxyz(ComovingMostBoundPosition, snapshot.GetComovingPosition(Elist[0].pid));
			copyHBTxyz(PhysicalMostBoundVelocity, snapshot.GetPhysicalVelocity(Elist[0].pid));
			break;
		  }
		}
	}
	ESnap.AverageKinematics(SpecificSelfPotentialEnergy, SpecificSelfKineticEnergy, SpecificAngularMomentum, Nbound, RefPos, RefVel);//only use CoM frame when unbinding and calculating Kinematics
	CountParticleTypes(snapshot);
#ifdef SAVE_BINDING_ENERGY
	Energies.resize(Nbound);
#pragma omp paralle for if(Nbound>100)
	for(HBTInt i=0;i<Nbound;i++)
	  Energies[i]=Elist[i].E;
#endif
}
void Subhalo_t::RecursiveUnbind(SubhaloList_t &Subhalos, const ParticleSnapshot_t &snap)
{
  bool is_orphan=(Nbound<=1);
  ParticleList_t particle_backup;
  if(is_orphan)	particle_backup=Particles;//orphans do not participate
  for(HBTInt i=0;i<NestedSubhalos.size();i++)
  {
    auto subid=NestedSubhalos[i];
	auto &subhalo=Subhalos[subid];
	subhalo.RecursiveUnbind(Subhalos, snap);
	Particles.insert(Particles.end(), subhalo.Particles.begin()+subhalo.Nbound, subhalo.Particles.end());
  }
  if(is_orphan)	Particles.swap(particle_backup);//restore to single particle for unbinding
  Unbind(snap);
  if(is_orphan)	Particles.swap(particle_backup);//set to extended list, to feed to its host
}

void Subhalo_t::TruncateSource()
{
  HBTInt Nsource;
  if(Nbound<=1) 
	Nsource=Nbound;
  else
	Nsource=Nbound*HBTConfig.SourceSubRelaxFactor;
  if(Nsource>Particles.size()) Nsource=Particles.size();
  Particles.resize(Nsource);
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
#ifdef INCLUSIVE_MASS
  #pragma omp parallel for schedule(dynamic,1) if(ParallelizeHaloes)
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	    Subhalos[subid].Unbind(*SnapshotPointer);
	    Subhalos[subid].TruncateSource();
  }
#else 
  HBTInt NumHalos=MemberTable.SubGroups.size();
  #pragma omp parallel for schedule(dynamic,1) if(ParallelizeHaloes)
    for(HBTInt haloid=0;haloid<NumHalos;haloid++)
    {
	  auto &subgroup=MemberTable.SubGroups[haloid];
	  if(subgroup.size()==0) continue;
	  //add new satellites to central's NestedSubhalos
	  auto &central=Subhalos[subgroup[0]];
	  auto &nests=central.NestedSubhalos;
	  auto old_membercount=nests.size();
	  auto &heads=MemberTable.SubGroupsOfHeads[haloid];
	  //update central member list (append other heads except itself)
	  nests.insert(nests.end(), heads.begin()+1, heads.end());
	  central.RecursiveUnbind(Subhalos, *SnapshotPointer);
	  nests.resize(old_membercount);//restore old satellite list
    }
  //unbind field subs  
  #pragma omp parallel
  {
      HBTInt NumField=MemberTable.SubGroups[-1].size();
    #pragma omp for schedule(dynamic,1) nowait
      for(HBTInt i=0; i<NumField;i++)
      {
	    HBTInt subid=MemberTable.SubGroups[-1][i];
	    Subhalos[subid].Unbind(*SnapshotPointer);
      }
    //unbind new-born subs
    HBTInt NumSubOld=MemberTable.AllMembers.size(), NumSub=Subhalos.size();
    #pragma omp for schedule(dynamic,1) 
      for(HBTInt i=NumSubOld;i<NumSub;i++)
      {
	    Subhalos[i].Unbind(*SnapshotPointer);
      }    
    #pragma omp for schedule(dynamic,1) 
      for(HBTInt i=0;i<NumSub;i++)
	    Subhalos[i].TruncateSource();
  }
#endif
}

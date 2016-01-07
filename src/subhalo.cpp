#include <iostream>
#include <new>
#include <algorithm>
#include <omp.h>

#include "datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"
#include "gravity_tree.h"

void MemberShipTable_t::SubIdToTrackId(const SubhaloList_t& Subhalos)
{
  /*not necessary currently*/
   for(HBTInt i=0;i<AllMembers.size();i++)
	AllMembers[i]=Subhalos[AllMembers[i]].TrackId;
}
void MemberShipTable_t::TrackIdToSubId(SubhaloList_t& Subhalos)
{
cout<<"Warning: TrackIdToSubId ToBe fully Implemented!\n";
// exit(1);
}

void SubhaloSnapshot_t::ParticleIdToIndex(const ParticleSnapshot_t& snapshot)
{//also bind to snapshot
#pragma omp single
  {
	SetEpoch(snapshot);
	SnapshotPointer=&snapshot;
  }
#pragma omp for
	for(HBTInt subid=0;subid<Subhalos.size();subid++)
	{
	  Subhalo_t::ParticleList_t & Particles=Subhalos[subid].Particles;
	  HBTInt nP=Particles.size();
	  for(HBTInt pid=0;pid<nP;pid++)
		Particles[pid]=snapshot.GetParticleIndex(Particles[pid]);
	}
}
void SubhaloSnapshot_t::ParticleIndexToId()
{
#pragma omp for
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	Subhalo_t::ParticleList_t & Particles=Subhalos[subid].Particles;
	HBTInt nP=Particles.size();
	for(HBTInt pid=0;pid<nP;pid++)
	  Particles[pid]=SnapshotPointer->GetParticleId(Particles[pid]);
  }
#pragma omp single
  SnapshotPointer=nullptr;
}

void SubhaloSnapshot_t::AverageCoordinates()
{
#pragma omp for
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	int coresize=GetCoreSize(Subhalos[subid].Nbound);
	SnapshotPointer->AveragePosition(Subhalos[subid].ComovingPosition, Subhalos[subid].Particles.data(), coresize);
	SnapshotPointer->AverageVelocity(Subhalos[subid].PhysicalVelocity, Subhalos[subid].Particles.data(), Subhalos[subid].Nbound);
  }
}

void Subhalo_t::CalculateProfileProperties(const ParticleSnapshot_t &part_snap)
{
  /* to calculate the following density-profile related properties
   *  currently only support const particle masses.
   * 
  HBTReal RmaxComoving;
  HBTReal VmaxPhysical;
  HBTReal LastMaxVmax;
  HBTInt SnapshotIndexOfLastMaxVmax; //the snapshot when it has the maximum Vmax, only considering past snapshots.
  
  HBTReal R2SigmaComoving;
  HBTReal RHalfComoving;
  
  HBTReal R200CritComoving;
  HBTReal R200MeanComoving;
  HBTReal RVirComoving;
  HBTReal M200Crit;
  HBTReal M200Mean;
  HBTReal MVir;
   */
  HBTReal PartMass=part_snap.GetMass(0);
  HBTReal VelocityUnit=sqrt(PhysicalConst::G*PartMass/part_snap.ScaleFactor);
  
  const HBTxyz &cen=part_snap.GetComovingPosition(Particles[0]); //most-bound particle as center.
  
  vector <HBTReal> r(Nbound), v(Nbound);
  for(HBTInt i=0;i<Nbound;i++)
	r[i]=PeriodicDistance(cen, part_snap.GetComovingPosition(Particles[i]));
  sort(r.begin(), r.end());
  for(HBTInt i=0;i<Nbound;i++)
  {
	  if(r[i]<HBTConfig.SofteningHalo) r[i]=HBTConfig.SofteningHalo; //resolution
	  v[i]=sqrt((HBTReal)(i+1)/r[i]);
  }
  HBTInt imax=max_element(v.begin(), v.end())-v.begin();
  RmaxComoving=r[imax];
  VmaxPhysical=v[imax]*VelocityUnit;
  RHalfComoving=r[Nbound/2];
  R2SigmaComoving=r[(HBTInt)(Nbound*0.955)];
  
  HBTReal virialF_tophat, virialF_b200, virialF_c200;
  part_snap.HaloVirialFactors(virialF_tophat, virialF_b200, virialF_c200);
  part_snap.SphericalOverdensitySize(MVir, RVirComoving, virialF_tophat, r);
  part_snap.SphericalOverdensitySize(M200Crit, R200CritComoving, virialF_c200, r);
  part_snap.SphericalOverdensitySize(M200Mean, R200MeanComoving, virialF_b200, r);
  
  if(VmaxPhysical>=LastMaxVmaxPhysical)
  {
	SnapshotIndexOfLastMaxVmax=part_snap.GetSnapshotIndex();
	LastMaxVmaxPhysical=VmaxPhysical;
  }

#ifdef ENABLE_EXPERIMENTAL_PROPERTIES
  /*the spin parameters are kind of ambiguous. do not provide*/
  for(int i=0;i<3;i++)
  {
	SpinPeebles[i]=SpecificAngularMomentum[i]*
	  sqrt(fabs(SpecificSelfPotentialEnergy+SpecificSelfKineticEnergy))/PhysicalConst::G/Nbound/PartMass;
	SpinBullock[i]=SpecificAngularMomentum[i]/sqrt(2.*PhysicalConst::G*PartMass*R2SigmaComoving);
  }
#endif
}
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
	Cosmology=snapshot.Cosmology;
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
#pragma omp parallel for
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
	Subhalo_t::ParticleList_t & Particles=Subhalos[subid].Particles;
	HBTInt nP=Particles.size();
	for(HBTInt pid=0;pid<nP;pid++)
	  Particles[pid]=SnapshotPointer->GetParticleId(Particles[pid]);
  }
  SnapshotPointer=nullptr;
}

void SubhaloSnapshot_t::AverageCoordinates()
{
#pragma omp for
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
// 	int coresize=GetCoreSize(Subhalos[subid].Nbound);
	copyHBTxyz(Subhalos[subid].ComovingMostBoundPosition, SnapshotPointer->GetComovingPosition(Subhalos[subid].Particles[0]));
	copyHBTxyz(Subhalos[subid].PhysicalMostBoundVelocity, SnapshotPointer->GetPhysicalVelocity(Subhalos[subid].Particles[0]));
	SnapshotPointer->AveragePosition(Subhalos[subid].ComovingAveragePosition, Subhalos[subid].Particles.data(), Subhalos[subid].Nbound);
	SnapshotPointer->AverageVelocity(Subhalos[subid].PhysicalAverageVelocity, Subhalos[subid].Particles.data(), Subhalos[subid].Nbound);
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
  if(Nbound<=1)
  {
	RmaxComoving=0.;
	VmaxPhysical=0.;
	R2SigmaComoving=0.;
	RHalfComoving=0.;
	R200CritComoving=0.;
	R200MeanComoving=0.;
	RVirComoving=0.;
	M200Crit=0.;
	M200Mean=0.;
	MVir=0.;
	for(int i=0;i<3;i++)
	{
	  SpinPeebles[i]=0.;
	  SpinBullock[i]=0.;
	}
	return;
  }
  HBTReal VelocityUnit=PhysicalConst::G/part_snap.Cosmology.ScaleFactor;
  
  const HBTxyz &cen=ComovingMostBoundPosition; //most-bound particle as center.
  
  vector <HBTReal> r(Nbound), v(Nbound);
  vector <double> m_in(Nbound);
  #pragma omp parallel if(Nbound>100)
  {
  #pragma omp for
  for(HBTInt i=0;i<Nbound;i++)
  {
	r[i]=PeriodicDistance(cen, part_snap.GetComovingPosition(Particles[i]));
	m_in[i]=part_snap.GetMass(Particles[i]);
  }
  #pragma omp single
  {
	sort(r.begin(), r.end());
	partial_sum(m_in.begin(), m_in.end(), m_in.begin());
  }
  #pragma omp for
  for(HBTInt i=0;i<Nbound;i++)
  {
	  if(r[i]<HBTConfig.SofteningHalo) r[i]=HBTConfig.SofteningHalo; //resolution
	  v[i]=m_in[i]/r[i];//v**2
  }
  }
  HBTInt imax=max_element(v.begin(), v.end())-v.begin();
  RmaxComoving=r[imax];
  VmaxPhysical=sqrt(v[imax]*VelocityUnit);
  RHalfComoving=r[Nbound/2];
  R2SigmaComoving=r[(HBTInt)(Nbound*0.955)];
  
  HBTReal virialF_tophat, virialF_b200, virialF_c200;
  part_snap.HaloVirialFactors(virialF_tophat, virialF_b200, virialF_c200);
  part_snap.SphericalOverdensitySize(MVir, RVirComoving, virialF_tophat, r, m_in);
  part_snap.SphericalOverdensitySize(M200Crit, R200CritComoving, virialF_c200, r, m_in);
  part_snap.SphericalOverdensitySize(M200Mean, R200MeanComoving, virialF_b200, r, m_in);
  
  if(VmaxPhysical>=LastMaxVmaxPhysical)
  {
	SnapshotIndexOfLastMaxVmax=part_snap.GetSnapshotIndex();
	LastMaxVmaxPhysical=VmaxPhysical;
  }

  HBTReal Mbound=GetMass();
  /*the spin parameters are kind of ambiguous. do not provide*/
  for(int i=0;i<3;i++)
  {
	SpinPeebles[i]=SpecificAngularMomentum[i]*
	  sqrt(fabs(SpecificSelfPotentialEnergy+SpecificSelfKineticEnergy))/PhysicalConst::G/Mbound;
	SpinBullock[i]=SpecificAngularMomentum[i]/sqrt(2.*PhysicalConst::G*Mbound*R2SigmaComoving);
  }
}

void Subhalo_t::CalculateShape(const ParticleSnapshot_t& part_snap)
{// (sum m_i r_i r_i/sum m_i)
  if(Nbound<=1)
  {
	#ifdef HAS_GSL
	for(int i=0;i<3;i++)
	  for(int j=0;j<3;j++)
	  {
		InertialEigenVector[i][j]=0.;
		InertialEigenVectorWeighted[i][j]=0.;
	  }
	#endif
	for(auto && I: InertialTensor) I=0.;
	for(auto && I: InertialTensorWeighted) I=0.;
  }
  const HBTxyz &cen=ComovingMostBoundPosition; //most-bound particle as center.
  
  double Ixx=0,Iyy=0, Izz=0, Ixy=0, Ixz=0, Iyz=0;
  double Ixxw=0,Iyyw=0, Izzw=0, Ixyw=0, Ixzw=0, Iyzw=0;
  #pragma omp parallel for reduction(+:Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Ixxw,Iyyw,Izzw,Ixyw,Ixzw,Iyzw) if(Nbound>100)
  for(HBTInt i=1;i<Nbound;i++)
  {
	  HBTReal m=part_snap.GetMass(Particles[i]);
	  const HBTxyz & pos=part_snap.GetComovingPosition(Particles[i]);
	  HBTReal dx=pos[0]-cen[0];
	  HBTReal dy=pos[1]-cen[1];
	  HBTReal dz=pos[2]-cen[2];
	  if(HBTConfig.PeriodicBoundaryOn)
	  {
		dx=NEAREST(dx);
		dy=NEAREST(dy);
		dz=NEAREST(dz);
	  }
	  HBTReal dx2=dx*dx;
	  HBTReal dy2=dy*dy;
	  HBTReal dz2=dz*dz;
	  Ixx+=dx2*m;
	  Iyy+=dy2*m;
	  Izz+=dz2*m;
	  Ixy+=dx*dy*m;
	  Ixz+=dx*dz*m;
	  Iyz+=dy*dz*m;
	  
	  HBTReal dr2=dx2+dy2+dz2;
	  dr2/=m;//for mass weighting
	  Ixxw+=dx2/dr2;
	  Iyyw+=dy2/dr2;
	  Izzw+=dz2/dr2;
	  Ixyw+=dx*dy/dr2;
	  Ixzw+=dx*dz/dr2;
	  Iyzw+=dy*dz/dr2;
  }
  InertialTensor[0]=Ixx; InertialTensor[1]=Ixy; InertialTensor[2]=Ixz; InertialTensor[3]=Iyy; InertialTensor[4]=Iyz; InertialTensor[5]=Izz;
  InertialTensorWeighted[0]=Ixxw; InertialTensorWeighted[1]=Ixyw; InertialTensorWeighted[2]=Ixzw; InertialTensorWeighted[3]=Iyyw; InertialTensorWeighted[4]=Iyzw; InertialTensorWeighted[5]=Izzw;
  HBTReal Mbound=GetMass();
  for(auto && I: InertialTensor) I/=Mbound;
  for(auto && I: InertialTensorWeighted) I/=Mbound;
#ifdef HAS_GSL  
  EigenAxis(Ixx, Ixy, Ixz, Iyy, Iyz, Izz, InertialEigenVector);
  EigenAxis(Ixxw, Ixyw, Ixzw, Iyyw, Iyzw, Izzw, InertialEigenVectorWeighted);
#endif
}

struct CompareParticleType_t
{
  const ParticleSnapshot_t & snap;
  CompareParticleType_t(const ParticleSnapshot_t & snap): snap(snap)
  {
  }
  bool operator () (const HBTInt & i, const HBTInt & j)
  {
	return snap.GetParticleType(i)<snap.GetParticleType(j);
  }
};

void Subhalo_t::CountParticleTypes(const ParticleSnapshot_t& part_snap)
{
//   CompareParticleType_t comparator(part_snap);
//   stable_sort(Particles.begin(), Particles.begin()+Nbound, comparator);//only sort the bound part
  for(int itype=0;itype<TypeMax;itype++)
  {
	NboundType[itype]=0;
	MboundType[itype]=0.;
  }
  if(Nbound>100)//parallelize
  {
	#pragma omp parallel
	{
	  vector <HBTInt> nboundtype(TypeMax, 0);
	  vector <float> mboundtype(TypeMax, 0.);
	  #pragma omp for
	  for(HBTInt i=0;i<Nbound;i++)
	  {
		auto &p=Particles[i];
		int itype=part_snap.GetParticleType(p);
		nboundtype[itype]++;
		mboundtype[itype]+=part_snap.GetMass(p);
	  }
	  #pragma omp critical
	  for(int i=0;i<TypeMax;i++)
	  {
		NboundType[i]+=nboundtype[i];
		MboundType[i]+=mboundtype[i];
	  }
	}
  }
  else
  {
	auto end=Particles.begin()+Nbound;
	for(auto it=Particles.begin();it!=end;++it)
	{
	  auto &p=*it;
	  int itype=part_snap.GetParticleType(p);
	  NboundType[itype]++;
	  MboundType[itype]+=part_snap.GetMass(p);
	}
  }
}


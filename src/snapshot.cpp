using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <chrono>

#include "snapshot.h"
#include "mymath.h"

class ParticleKeyList_t: public KeyList_t <HBTInt, HBTInt>
{
  typedef HBTInt Index_t;
  typedef HBTInt Key_t;
  const ParticleSnapshot_t &Snap;
public:
  ParticleKeyList_t(ParticleSnapshot_t &snap): Snap(snap) {};
  Index_t size() const
  {
        return Snap.size();
  }
  Key_t GetKey(Index_t i) const
  {
        return Snap.GetMemberId(i);
  }
  Index_t GetIndex(Index_t i) const
  {
        return i;
  }
};

void ParticleSnapshot_t::FillParticleHash()
{
  ParticleKeyList_t Ids(*this);
  ParticleHash->Fill(Ids);
}

void ParticleSnapshot_t::ClearParticleHash()
{
  ParticleHash->Clear();
}

void ParticleSnapshot_t::Clear()
/*reset to empty*/
{
  vector<Particle_t>().swap(Particles);
  ClearParticleHash();//even if you don't do this, the destructor will still clean up the memory.
}

double ParticleSnapshot_t::AveragePosition(HBTxyz& CoM, const ParticleIndex_t PartIndex[], const ParticleIndex_t NumPart) const
/*mass weighted average position*/
{
	ParticleIndex_t i,j;
	double sx[3],origin[3],msum;
	
	if(0==NumPart) return 0.;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoM, GetComovingPosition(PartIndex[0]));
	  return GetParticleMass(PartIndex[0]);
	}
	
	sx[0]=sx[1]=sx[2]=0.;
	msum=0.;
	if(HBTConfig.PeriodicBoundaryOn)
	  for(j=0;j<3;j++)
		origin[j]=GetComovingPosition(PartIndex[0])[j];
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=GetParticleMass(PartIndex[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
	  if(HBTConfig.PeriodicBoundaryOn)
		  sx[j]+=NEAREST(GetComovingPosition(PartIndex[i])[j]-origin[j])*m;
	  else
		  sx[j]+=GetComovingPosition(PartIndex[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	{
		sx[j]/=msum;
		if(HBTConfig.PeriodicBoundaryOn) sx[j]+=origin[j];
		CoM[j]=sx[j];
	}
	return msum;
}
double ParticleSnapshot_t::AverageVelocity(HBTxyz& CoV, const ParticleIndex_t PartIndex[], const ParticleIndex_t NumPart) const
/*mass weighted average velocity*/
{
	ParticleIndex_t i,j;
	double sv[3],msum;
	
	if(0==NumPart) return 0.;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoV, GetPhysicalVelocity(PartIndex[0]));
	  return GetParticleMass(PartIndex[0]);
	}
	
	sv[0]=sv[1]=sv[2]=0.;
	msum=0.;
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=GetParticleMass(PartIndex[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
		sv[j]+=GetPhysicalVelocity(PartIndex[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	  CoV[j]=sv[j]/msum;
	return msum;
}

//TODO: detach these SO functions from snapshot. pass cosmology as parameter
void Cosmology_t::SphericalOverdensitySize(float& Mvir, float& Rvir, HBTReal VirialFactor, const vector< HBTReal >& RSorted) const
/*
 * find SphericalOverdensitySize from a given list of sorted particle distances with default ParticleMass in Comology.
 * all distances comoving.
 * 
 * Brute-force method to scan from most-distant particle.
 * 
 */
{
  assert(ParticleMass>0);
  HBTInt i, np=RSorted.size();
  HBTReal RhoVirial=VirialFactor*Hz*Hz/2.0/PhysicalConst::G/ParticleMass*ScaleFactor*ScaleFactor*ScaleFactor;
  for(i=np;i>0;i--)
  {
	HBTReal r=RSorted[i-1];
	if(i>RhoVirial*r*r*r) break;
  }
  Mvir=i*ParticleMass;
  Rvir=pow(i/RhoVirial, 1.0/3);//comoving
}

void Cosmology_t::SphericalOverdensitySize(float& Mvir, float& Rvir, HBTReal VirialFactor, const vector <RadVelMass_t> &prof) const
/*
 * find SphericalOverdensitySize from a given list of sorted particle distances.
 * all distances comoving.
 * 
 * Brute-force method to scan from most-distant particle.
 * 
 */
{ 
  HBTReal RhoVirial=VirialFactor*Hz*Hz/2.0/PhysicalConst::G*ScaleFactor*ScaleFactor*ScaleFactor;
  
  for(auto p=prof.rbegin();p!=prof.rend();++p)
  {
	auto r=p->r;
	auto m=p->m;
	if(m>RhoVirial*r*r*r) 
	{
	  Mvir=m;
	  Rvir=pow(m/RhoVirial, 1.0/3);//comoving
	  return;
	}
  }
}

void Cosmology_t::SphericalOverdensitySize2(float& Mvir, float& Rvir, HBTReal VirialFactor, const vector< HBTReal >& RSorted) const
/*
 * find SphericalOverdensitySize from a given list of sorted particle distances with default ParticleMass in Comology.
 * all distances comoving.
 * 
 * Iterative method to guess the radius.
 * 
 */
{
  assert(ParticleMass>0);
  HBTReal tol=1e-5;
  HBTInt i,ndiv, np=RSorted.size();
  HBTReal rvir,rdiv;
  
  ndiv=np;//guess mass
  rdiv=RSorted[ndiv-1];
  rvir=pow(2.0*PhysicalConst::G*ndiv*ParticleMass/VirialFactor/Hz/Hz,1.0/3)/ScaleFactor;//guess radius
  while(true)
  {
	if(rdiv>rvir)//reduce mass guess
	{
	  for(i=ndiv-1;i>=0;i--)
	  {
		if(RSorted[i]<rvir) break;
	  }
	  ndiv=i+1;
	}
	else if(rdiv<rvir) //increase mass guess
	{
	  for(i=ndiv;i<np;i++)
	  {
		if(RSorted[i]>=rvir) break;
	  }
	  ndiv=i;
	}
	
	rdiv=rvir;
	rvir=pow(2.0*PhysicalConst::G*ndiv*ParticleMass/VirialFactor/Hz/Hz,1.0/3)/ScaleFactor;//recalibrate radius
	
	if(0==ndiv||np==ndiv||fabs((rvir-rdiv)/rvir)<tol) break; //converged
  }
  
  Rvir=rvir;
  Mvir=ndiv*ParticleMass;
}

void Cosmology_t::HaloVirialFactors(HBTReal &virialF_tophat, HBTReal &virialF_b200, HBTReal &virialF_c200) const
//obtain <Rho_vir>/Rho_cri
{
	HBTReal Hratio,x;
	Hratio=Hz/PhysicalConst::H0;
	x=OmegaZ-1;
	virialF_tophat=18.0*3.1416*3.1416+82.0*x-39.0*x*x;//<Rho_vir>/Rho_cri
	virialF_c200=200.;
	virialF_b200=200.*OmegaZ;//virialF w.r.t contemporary critical density 
}
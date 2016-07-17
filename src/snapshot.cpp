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
#include <mpi.h>

ostream& operator << (ostream& o, Particle_t &p)
{
   o << "[" << p.Id << ", " <<p.Mass<<", " << p.ComovingPosition << ", " << p.PhysicalVelocity << "]";
   return o;
};
/*
void ParticleSnapshot_t::FillParticleHash()
{
  cout<<"Filling Hash Table...\n";
  auto begin = chrono::high_resolution_clock::now();
  ParticleHash.rehash(NumberOfParticles);
  ParticleHash.reserve(NumberOfParticles);
  for(HBTInt i=0;i<NumberOfParticles;i++)
	ParticleHash[ParticleId[i]]=i;
  auto end = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::milliseconds>(end - begin);
  cout << "inserts: " << elapsed.count() << endl;
//   cout<<ParticleHash.bucket_count()<<" buckets used; load factor "<<ParticleHash.load_factor()<<endl;
}
void ParticleSnapshot_t::ClearParticleHash()
{
  ParticleHash.clear();
}
*/
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
	return Snap.GetId(i);
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
//   ParticleHash->GetKeyMinMax(IdMin, IdMax);
}
void ParticleSnapshot_t::ClearParticleHash()
{
  ParticleHash->Clear();
}

void ParticleSnapshot_t::AveragePosition(HBTxyz& CoM, const HBTInt Particles[], HBTInt NumPart) const
/*mass weighted average position*/
{
	HBTInt i,j;
	double sx[3],origin[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoM, GetComovingPosition(Particles[0]));
	  return;
	}
	
	sx[0]=sx[1]=sx[2]=0.;
	msum=0.;
	if(HBTConfig.PeriodicBoundaryOn)
	  for(j=0;j<3;j++)
		origin[j]=GetComovingPosition(Particles[0])[j];
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=GetMass(Particles[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
	  if(HBTConfig.PeriodicBoundaryOn)
		  sx[j]+=NEAREST(GetComovingPosition(Particles[i])[j]-origin[j])*m;
	  else
		  sx[j]+=GetComovingPosition(Particles[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	{
		sx[j]/=msum;
		if(HBTConfig.PeriodicBoundaryOn) sx[j]+=origin[j];
		CoM[j]=sx[j];
	}
}
void ParticleSnapshot_t::AverageVelocity(HBTxyz& CoV, const HBTInt Particles[], HBTInt NumPart) const
/*mass weighted average velocity*/
{
	HBTInt i,j;
	double sv[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoV, GetPhysicalVelocity(Particles[0]));
	  return;
	}
	
	sv[0]=sv[1]=sv[2]=0.;
	msum=0.;
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=GetMass(Particles[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
		sv[j]+=GetPhysicalVelocity(Particles[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	  CoV[j]=sv[j]/msum;
}

void Particle_t::create_MPI_type(MPI_Datatype &dtype)
{
/*to create the struct data type for communication*/	
Particle_t &p=*this;
#define MaxNumAttr 10
MPI_Datatype oldtypes[MaxNumAttr];
MPI_Aint   offsets[MaxNumAttr], origin,extent;
int blockcounts[MaxNumAttr];

MPI_Get_address(&p,&origin);
MPI_Get_address((&p)+1,&extent);//to get the extent of s
extent-=origin;
  
  int i=0;
  #define RegisterAttr(x, type, count) {MPI_Get_address(&(p.x), offsets+i); offsets[i]-=origin; oldtypes[i]=type; blockcounts[i]=count; i++;}
  RegisterAttr(Id, MPI_HBT_INT, 1)
  RegisterAttr(ComovingPosition, MPI_HBT_REAL, 3)
  RegisterAttr(PhysicalVelocity, MPI_HBT_REAL, 3)
#ifndef DM_ONLY
  RegisterAttr(Mass, MPI_HBT_REAL, 1)
#ifdef UNBIND_WITH_THERMAL_ENERGY
  RegisterAttr(InternalEnergy, MPI_HBT_REAL, 1)
#endif
  RegisterAttr(Type, MPI_INT, 1)
#endif
  #undef RegisterAttr
  assert(i<=MaxNumAttr);
  
  MPI_Type_create_struct(i,blockcounts,offsets,oldtypes, &dtype);
  MPI_Type_create_resized(dtype,(MPI_Aint)0, extent, &dtype);
  MPI_Type_commit(&dtype);
  #undef MaxNumAttr
}

double AveragePosition(HBTxyz& CoM, const Particle_t Particles[], HBTInt NumPart)
/*mass weighted average position*/
{
	HBTInt i,j;
	double sx[3],origin[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoM, Particles[0].ComovingPosition);
	  return;
	}
	
	sx[0]=sx[1]=sx[2]=0.;
	msum=0.;
	if(HBTConfig.PeriodicBoundaryOn)
	  for(j=0;j<3;j++)
		origin[j]=Particles[0].ComovingPosition[j];
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=Particles[i].Mass;
	  msum+=m;
	  for(j=0;j<3;j++)
	  if(HBTConfig.PeriodicBoundaryOn)
		  sx[j]+=NEAREST(Particles[i].ComovingPosition[j]-origin[j])*m;
	  else
		  sx[j]+=Particles[i].ComovingPosition[j]*m;
	}
	
	for(j=0;j<3;j++)
	{
		sx[j]/=msum;
		if(HBTConfig.PeriodicBoundaryOn) sx[j]+=origin[j];
		CoM[j]=sx[j];
	}
	return msum;
}
double AverageVelocity(HBTxyz& CoV, const Particle_t Particles[], HBTInt NumPart)
/*mass weighted average velocity*/
{
	HBTInt i,j;
	double sv[3],msum;
	
	if(0==NumPart) return;
	if(1==NumPart) 
	{
	  copyHBTxyz(CoV, Particles[0].PhysicalVelocity);
	  return;
	}
	
	sv[0]=sv[1]=sv[2]=0.;
	msum=0.;
	
	for(i=0;i<NumPart;i++)
	{
	  HBTReal m=Particles[i].Mass;
	  msum+=m;
	  for(j=0;j<3;j++)
		sv[j]+=Particles[i].PhysicalVelocity[j]*m;
	}
	
	for(j=0;j<3;j++)
	  CoV[j]=sv[j]/msum;
	
	return msum;
}

void Snapshot_t::SphericalOverdensitySize(float& Mvir, float& Rvir, HBTReal VirialFactor, const vector< HBTReal >& RSorted, HBTReal ParticleMass) const
/*
 * find SphericalOverdensitySize from a given list of sorted particle distances.
 * all distances comoving.
 * 
 * Brute-force method to scan from most-distant particle.
 * 
 * currently only support constant particle mass.
 */
{  
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

void Snapshot_t::SphericalOverdensitySize2(float& Mvir, float& Rvir, HBTReal VirialFactor, const vector< HBTReal >& RSorted, HBTReal ParticleMass) const
/*
 * find SphericalOverdensitySize from a given list of sorted particle distances.
 * all distances comoving.
 * 
 * Iterative method to guess the radius.
 * 
 * currently only support constant particle mass.
 */
{
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

void Snapshot_t::HaloVirialFactors(HBTReal &virialF_tophat, HBTReal &virialF_b200, HBTReal &virialF_c200) const
{
	HBTReal Hratio,x,OmegaZ;
	Hratio=Hz/PhysicalConst::H0;
	OmegaZ=OmegaM0/(ScaleFactor*ScaleFactor*ScaleFactor)/Hratio/Hratio;
	x=OmegaZ-1;
	virialF_tophat=18.0*3.1416*3.1416+82.0*x-39.0*x*x;//<Rho_vir>/Rho_cri
	virialF_c200=200.;
	virialF_b200=200.*OmegaZ;//virialF w.r.t contemporary critical density 
}

void ParticleSnapshot_t::Clear()
/*reset to empty*/
{
  #define RESET(x, T) {vector <T>().swap(x);}
  RESET(Particles, Particle_t);
  #undef RESET 
  ClearParticleHash();//even if you don't do this, the destructor will still clean up the memory.
  //   cout<<NumberOfParticles<<" particles cleared from snapshot "<<SnapshotIndex<<endl;
  NumberOfParticlesOnAllNodes=0;
}

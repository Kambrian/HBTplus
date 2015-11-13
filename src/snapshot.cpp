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

void ParticleSnapshot_t::FillParticleHash()
{
  if(HBTConfig.ParticleIdNeedHash)
	ParticleHash=&MappedHash;
  else
	ParticleHash=&FlatHash;
  cout<<"Filling Hash Table...\n";
  auto begin = chrono::high_resolution_clock::now();
  ParticleHash->Fill(ParticleId.data(), NumberOfParticles);
  auto end = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<chrono::duration<double>>(end - begin);
  cout << "HashTableFilled in: " << elapsed.count() <<"seconds"<< endl;
}
void ParticleSnapshot_t::ClearParticleHash()
{
  ParticleHash->Clear();
}

void ParticleSnapshot_t::AveragePosition(HBTxyz& CoM, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const
/*mass weighted average position*/
{
	ParticleIndex_t i,j;
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
	  HBTReal m=GetParticleMass(Particles[i]);
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
void ParticleSnapshot_t::AverageVelocity(HBTxyz& CoV, const ParticleIndex_t Particles[], const ParticleIndex_t NumPart) const
/*mass weighted average velocity*/
{
	ParticleIndex_t i,j;
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
	  HBTReal m=GetParticleMass(Particles[i]);
	  msum+=m;
	  for(j=0;j<3;j++)
		sv[j]+=GetPhysicalVelocity(Particles[i])[j]*m;
	}
	
	for(j=0;j<3;j++)
	  CoV[j]=sv[j]/msum;
}

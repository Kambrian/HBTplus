#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <glob.h>
#include <climits>
#include <algorithm>
#include <chrono>

#include "mymath.h"
#include "halo.h"

void HaloSnapshot_t::ParticleIdToIndex(const ParticleSnapshot_t& snapshot)
{
  chrono::high_resolution_clock::time_point time_begin, time_end;
#pragma omp master
  {
  SetEpoch(snapshot);
  ParticleSnapshot=&snapshot;
  time_begin = chrono::high_resolution_clock::now();
  }
#pragma omp for //maybe try collapse(2)? need to remove intermediate variables to enable this.
  for(HBTInt haloid=0;haloid<Halos.size();haloid++)
  {
	Halo_t::ParticleList_t & Particles=Halos[haloid].Particles;
	HBTInt np=Particles.size();
	for(HBTInt pid=0;pid<np;pid++)
	  Particles[pid]=snapshot.GetParticleIndex(Particles[pid]);//this should be safe since its a const func
  }
  #pragma omp master
  {
	time_end = chrono::high_resolution_clock::now();
	auto elapsed = chrono::duration_cast<chrono::duration<double>>(time_end - time_begin);
	cout << "Halo ParticleId lookup took " << elapsed.count() <<"seconds"<< endl;
  }
}

void HaloSnapshot_t::ParticleIndexToId()
{
#pragma omp for
  for(HBTInt haloid=0;haloid<Halos.size();haloid++)
  {
	Halo_t::ParticleList_t &Particles=Halos[haloid].Particles;
	HBTInt nP=Halos[haloid].Particles.size();
	for(HBTInt pid=0;pid<nP;pid++)
	  Particles[pid]=ParticleSnapshot->GetParticleId(Particles[pid]);
  }
#pragma omp single
  ParticleSnapshot=nullptr;
}

void HaloSnapshot_t::AverageCoordinates()
{
#pragma omp for
  for(HBTInt i=0;i<Halos.size();i++)
  {
	ParticleSnapshot->AveragePosition(Halos[i].ComovingAveragePosition, Halos[i].Particles.data(), Halos[i].Particles.size());
	ParticleSnapshot->AverageVelocity(Halos[i].PhysicalAverageVelocity, Halos[i].Particles.data(), Halos[i].Particles.size());
  }
}
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

#include "../snapshot.h"
#include "../mymath.h"
#include "gadget_io.h"

void ParticleSnapshot_t::Load(int snapshot_index, bool fill_particle_hash)
{ 
  SetSnapshotIndex(snapshot_index);
  
  if(HBTConfig.SnapshotFormat=="gadget")
  {
	GadgetReader_t(SnapshotId, Particles, Cosmology);
  }
  
  if(HBTConfig.SnapshotFormat=="mysnapshot")
  {//insert your snapshot reader here, and include relevant header in the header if necessary
  }
  
  if(fill_particle_hash)
	FillParticleHash();
}


#ifdef TEST_snapshot_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  HBTConfig.ParseConfigFile(argv[1]);
  ParticleSnapshot_t snapshot;
  snapshot.Load(HBTConfig.MaxSnapshotIndex, true);
  cout<<snapshot.GetNumberOfParticles()<<endl;
  cout<<snapshot.GetParticleId(10)<<endl;
  cout<<snapshot.GetComovingPosition(10)<<endl;
  cout<<snapshot.GetParticleMass(10)<<','<<snapshot.GetParticleMass(100)<<endl;
  cout<<snapshot.GetParticleIndex(snapshot.GetParticleId(10))<<endl;
  return 0;
}
#endif
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
#include "apostle_io.h"
#ifdef USE_JING_IO
#include "jing/jing_io.h"
#endif

void ParticleSnapshot_t::Load(int snapshot_index, bool fill_particle_hash)
{
  SetSnapshotIndex(snapshot_index);

  if(HBTConfig.SnapshotFormat=="gadget")
  {
	GadgetReader_t(SnapshotId, Particles, Cosmology);
  }
  else if(Apostle::IsIllustrisSnap(HBTConfig.SnapshotFormat)||
      Apostle::IsHelucidSnap(HBTConfig.SnapshotFormat)||
      Apostle::IsApostleSnap(HBTConfig.SnapshotFormat))
  {
	Apostle::ApostleReader_t().LoadSnapshot(SnapshotId, Particles, Cosmology);
  }
#ifdef USE_JING_IO
  else if(HBTConfig.SnapshotFormat=="jing")
  {
    JingReader_t reader(SnapshotId);
    reader.LoadSnapshot(Particles, Cosmology);
  }
#endif
  else if(HBTConfig.SnapshotFormat=="mysnapshot")
  {/*insert your snapshot reader here, and include relevant header in the header if necessary
	you need to fill up Particles vector, and set the cosmology, e.g.,

	LoadMySnapshot(SnapshotId, Particles, Cosmology);

	*/
  }
  else
	throw(runtime_error("unknown SnapshotFormat "+HBTConfig.SnapshotFormat));

#ifdef DM_ONLY
  assert(Cosmology.ParticleMass>0);
#endif

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
  cout<<"ParticleMass:"<<snapshot.Cosmology.ParticleMass<<endl;
  cout<<"Number Of Particles:"<<snapshot.GetNumberOfParticles()<<endl;
  cout<<"Id of 10th particle: "<<snapshot.GetParticleId(10)<<endl;
  cout<<"position: "<<snapshot.GetComovingPosition(10)<<endl;
  cout<<"Mass: "<<snapshot.GetParticleMass(10)<<", Mass of 100th particle:"<<snapshot.GetParticleMass(100)<<endl;
  cout<<"Index of 10th particle : "<<snapshot.GetParticleIndex(snapshot.GetParticleId(10))<<endl;
  return 0;
}
#endif

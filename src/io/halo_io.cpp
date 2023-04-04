#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <glob.h>
#include <climits>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstdio>

#include "../mymath.h"
#include "../halo.h"
#include "gadget_group_io.h"
#include "apostle_io.h"

void HaloSnapshot_t::Load(MpiWorker_t &world, int snapshot_index)
{
  SetSnapshotIndex(snapshot_index);

  string GroupFileFormat=HBTConfig.GroupFileFormat;

  if(GadgetGroup::IsGadgetGroup(GroupFileFormat))
	GadgetGroup::Load(world, SnapshotId, Halos);
  else if(IsApostleGroup(GroupFileFormat))
	ApostleReader_t().LoadGroups(world, SnapshotId, Halos);
  else if(GroupFileFormat=="my_group_format")
  {/*extend your own group reader here, input SnapshotId and output filled Halo list, e.g.:

	MyGroupReader(world, SnapshotId, Halos)

	*/
  }
  else
	throw(runtime_error("unknown GroupFileFormat "+GroupFileFormat));

  NumPartOfLargestHalo=0;
  TotNumberOfParticles=0;
  for(auto && h: Halos)
  {
    auto np=h.Particles.size();
    TotNumberOfParticles+=np;
    if(np>NumPartOfLargestHalo) NumPartOfLargestHalo=np;
  }

  HBTInt NumHalos=Halos.size(), NumHalosAll=0;
  MPI_Reduce(&NumHalos, &NumHalosAll, 1, MPI_HBT_INT, MPI_SUM, 0, world.Communicator);
  if(world.rank()==0)
    cout<<NumHalosAll<<" groups loaded at snapshot "<<snapshot_index<<"("<<SnapshotId<<")"<<endl;
}

void HaloSnapshot_t::LoadWithSnap(MpiWorker_t &world, int snapshot_index, ParticleSnapshot_t &partsnap)
{
  SetSnapshotIndex(snapshot_index);

  string GroupFileFormat=HBTConfig.GroupFileFormat;

  if(Gadget4Group::IsGadget4Group(GroupFileFormat))
	Gadget4Group::Load(world, SnapshotId, Halos, partsnap);
  else
	throw(runtime_error("unknown GroupFileFormat to LoadWithSnap: "+GroupFileFormat));

  NumPartOfLargestHalo=0;
  TotNumberOfParticles=0;
  for(auto && h: Halos)
  {
    auto np=h.Particles.size();
    TotNumberOfParticles+=np;
    if(np>NumPartOfLargestHalo) NumPartOfLargestHalo=np;
  }

  HBTInt NumHalos=Halos.size(), NumHalosAll=0;
  MPI_Reduce(&NumHalos, &NumHalosAll, 1, MPI_HBT_INT, MPI_SUM, 0, world.Communicator);
  if(world.rank()==0)
    cout<<NumHalosAll<<" groups loaded at snapshot "<<snapshot_index<<"("<<SnapshotId<<")"<<endl;
}

#ifdef TEST_halo_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MpiWorker_t world(MPI_COMM_WORLD);

  int snapshot_start, snapshot_end;
  if(0==world.rank())
	ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  HBTConfig.BroadCast(world, 0, snapshot_start, snapshot_end);

  HaloSnapshot_t halo;
//   ParticleSnapshot_t snap;
//   snap.Load(world, snapshot_start);
//   cout<<"snapshot loaded\n";
  halo.Load(world, snapshot_start);
//   halo.UpdateParticles(world, snap);
  if(halo.Halos.size()>1)
  {
	auto & h=halo.Halos[1];
	cout<<" Halo 1 from thread "<<world.rank()<<":"<<"id="<<h.HaloId<<","<<h.Particles.size()<<", "<<h.Particles[5]<<endl;
  }

  MPI_Finalize();
  return 0;
}
#endif

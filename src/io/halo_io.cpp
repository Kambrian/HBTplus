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
#include <cmath>

#include "../mymath.h"
#include "../halo.h"
#include "gadget_group_io.h"
#include "apostle_io.h"
#include "gadget4_io.h"

void HaloSnapshot_t::Load ( MpiWorker_t& world, int snapshot_index )
{//compatibility interface
  if(Gadget4Reader::IsGadget4Group(HBTConfig.GroupFileFormat))
  {
    throw(runtime_error("reading gadget4 group requires preloaded snapshot\n"));
    exit(1);
  }

  ParticleSnapshot_t partsnap;
  partsnap.SetSnapshotIndex(snapshot_index);
  Load(world, partsnap);
}

void HaloSnapshot_t::Load(MpiWorker_t &world, const ParticleSnapshot_t &partsnap)
{
  int snapshot_index=partsnap.GetSnapshotIndex();
  SetSnapshotIndex(snapshot_index);

  string GroupFileFormat=HBTConfig.GroupFileFormat;


  if(GadgetGroup::IsGadgetGroup(GroupFileFormat))
  {
    if(Gadget4Reader::IsGadget4Group(GroupFileFormat))
      Gadget4Reader::Gadget4Reader_t().LoadGroups(world, partsnap, Halos);//only this needs partsnap
    else
      GadgetGroup::Load(world, SnapshotId, Halos);
  }
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
 #pragma omp parallel for reduction(max:NumPartOfLargestHalo) reduction(+:TotNumberOfParticles)
  for(HBTInt i=0;i<Halos.size();i++)
  {
    HBTInt np=Halos[i].Particles.size();
    TotNumberOfParticles+=np;
    if(np>NumPartOfLargestHalo) NumPartOfLargestHalo=np;
  }

  HBTInt NumHalos=Halos.size(), NumHalosAll=0;
  MPI_Reduce(&NumHalos, &NumHalosAll, 1, MPI_HBT_INT, MPI_SUM, 0, world.Communicator);
  if(world.rank()==0)
    cout<<NumHalosAll<<" groups loaded at snapshot "<<snapshot_index<<"("<<SnapshotId<<")"<<endl;

  HBTInt MaxSizeAll=0;
  MPI_Reduce(&NumPartOfLargestHalo, &MaxSizeAll, 1, MPI_HBT_INT, MPI_MAX, 0, world.Communicator);
  if(world.rank()==0) cout<<"Max halo size="<<MaxSizeAll<<endl;
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

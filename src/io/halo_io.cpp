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

#include "../mymath.h"
#include "../halo.h"
#include "gadget_group_io.h"
#include "apostle_io.h"
#include "jing/jing_io.h"
#include "custom_io.h"

void HaloSnapshot_t::Load(int snapshot_index)
{
  SetSnapshotIndex(snapshot_index);
  
  string GroupFileFormat=HBTConfig.GroupFileFormat;
  
  if(GadgetGroup::IsGadgetGroup(GroupFileFormat))
	TotNumberOfParticles=GadgetGroup::Load(SnapshotId, Halos);
  else if(IsApostleGroup(GroupFileFormat))
	TotNumberOfParticles=ApostleReader_t().LoadGroups(SnapshotId, Halos);
  else if(JingGroup::IsJingGroup(HBTConfig.GroupFileFormat))
	TotNumberOfParticles=JingGroup::LoadGroup(SnapshotId, Halos);
  else if(GroupFileFormat=="my_group_format")
  {/*extend your own group reader here, fill Halos and return TotNumberOfParticles, e.g.:*/
	TotNumberOfParticles=MyGroupReader(SnapshotId, Halos);
  }
  else
	throw(runtime_error("unknown GroupFileFormat "+GroupFileFormat));
  
  NumPartOfLargestHalo=0;
  for(auto && h: Halos)
	if(h.Particles.size()>NumPartOfLargestHalo) NumPartOfLargestHalo=h.Particles.size();
}


#ifdef TEST_halo_io
#include "../config_parser.h"

int main(int argc, char **argv)
{
  HBTConfig.ParseConfigFile(argv[1]);
  HaloSnapshot_t halo;
  halo.Load(HBTConfig.MaxSnapshotIndex);
  cout<<halo.Halos.size()<<";"<<halo.TotNumberOfParticles<<endl;
  cout<<halo.Halos[10].Particles.size()<<endl;
  cout<<halo.Halos[10].Particles[0]<<endl;
/*  ParticleSnapshot_t snap;
  snap.Load(HBTConfig.MaxSnapshotIndex);
  halo.ParticleIdToIndex(snap);
  int ids[]={1,53};
  for(auto &hid: ids)
  {
	auto & h=halo.Halos[hid];
	auto pid=h.Particles[5];
	const HBTxyz &pos=snap.GetComovingPosition(pid);
	const HBTxyz &vel=snap.GetPhysicalVelocity(pid);
	cout<<" Halo id="<<hid<<","<<h.Particles.size()<<", ["<<snap.GetParticleId(pid)<<","<<snap.GetParticleMass(pid)<<",("<<pos[0]<<","<<pos[1]<<","<<pos[2]<<"), ("<<vel[0]<<","<<vel[1]<<","<<vel[2]<<") ]"<<endl;
  }
  */
  return 0;
}
#endif
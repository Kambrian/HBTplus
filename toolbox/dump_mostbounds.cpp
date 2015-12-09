using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "../src/snapshot.h"
#include "../src/subhalo.h"
#include "../src/mymath.h"

int main(int argc, char **argv)
{
  int snapshot_start, snapshot_end;
  ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  for(int isnap=snapshot_start;isnap<=snapshot_end;isnap++)
  {
  SubhaloSnapshot_t subsnap;
  subsnap.Load(isnap);
  ParticleSnapshot_t partsnap;
  partsnap.Load(isnap);
  subsnap.ParticleIdToIndex(partsnap);
  
  for(int subid=0;subid<2;subid++)
  {  
	stringstream filename;
	filename<<HBTConfig.SubhaloPath<<"/postproc/Particles_"<<isnap<<"."<<subid;
	ofstream outfile(filename.str().c_str(), fstream::out);
	if(subsnap.Subhalos[subid].Particles.size()>10)
	  for(int i=0;i<10;i++)
	  {
		HBTInt pid=subsnap.Subhalos[subid].Particles[i];
		for(auto &&x: partsnap.GetComovingPosition(pid))
		  outfile<<x<<" ";
		outfile<<endl;
	  }
  }
  }
  
  return 0;
}

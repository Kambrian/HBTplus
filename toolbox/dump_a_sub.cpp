using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "../simulation_io/snapshot.h"
#include "../simulation_io/subhalo.h"
#include "../mymath.h"

int main(int argc, char **argv)
{
  int isnap=32, subid=22;
  HBTConfig.ParseConfigFile(argv[1]);
  SubhaloSnapshot_t subsnap;
  subsnap.Load(isnap, true);

  FILE *fp;
  stringstream filename;
  filename<<HBTConfig.SubhaloPath<<"/postproc/Subhalo_"<<isnap<<"."<<subid;
  myfopen(fp, filename.str().c_str(), "w");
  cout<<subsnap.Subhalos[subid].Particles.size()<<endl;
  fwrite(subsnap.Subhalos[subid].Particles.data(), sizeof(HBTInt), subsnap.Subhalos[subid].Particles.size(), fp);
  fclose(fp);
  
  return 0;
}

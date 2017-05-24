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
  HBTConfig.ParseConfigFile(argv[1]);
  int isnap=atoi(argv[2]), subid=atoi(argv[3]);
  SubhaloSnapshot_t subsnap;
  subsnap.Load(isnap);

  cout<<subsnap.Subhalos[subid].Particles.size()<<endl;
  
  stringstream filename;
  filename<<HBTConfig.SubhaloPath<<"/postproc";
  mkdir(filename.str().c_str(), 0755);
  filename<<"/Subhalo_"<<isnap<<"."<<subid;
  ofstream outfile;
  outfile.open(filename.str(), fstream::out|fstream::app);
  for(auto &&p: subsnap.Subhalos[subid].Particles)
    outfile<<p<<endl;
//   FILE *fp;
//   myfopen(fp, filename.str().c_str(), "w");
//   fwrite(subsnap.Subhalos[subid].Particles.data(), sizeof(HBTInt), subsnap.Subhalos[subid].Particles.size(), fp);
//   fclose(fp);
  
  return 0;
}

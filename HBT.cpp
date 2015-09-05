using namespace std;
#include <iostream>
// #include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "simulation_io/snapshot.h"
#include "simulation_io/halo.h"
#include "simulation_io/subhalo.h"
#include "mymath.h"

int main(int argc, char **argv)
{
  int snapshot_start, snapshot_end;
  ParseHBTParams(argc, argv, HBTConfig, snapshot_start, snapshot_end);
  SubHaloSnapshot_t subsnap;
  
  subsnap.Load(HBTConfig, snapshot_start-1);
  
  for(int isnap=snapshot_start;isnap<=snapshot_end;isnap++)
  {
	Snapshot_t snapshot;
	HaloSnapshot_t halosnap;
	snapshot.Load(HBTConfig, isnap);
	halosnap.Load(HBTConfig, isnap);
	
	subsnap.rehost(halosnap);
	subsnap.filter();
	
	subsnap.save();
// 	subsnap_old=subsnap;
  }
  
  return 0;
}

void ParseHBTParams(int argc, char **argv, Parameter_t &config, int &snapshot_start, int &snapshot_end)
{
  if(argc<2)
  {
	cerr<<"Usage: "<<argv[0]<<" [param_file] <snapshot_start> <snapshot_end>\n";
	exit(1);
  }
  config.ParseConfigFile(argv[1]);
  if(argc>2)
	snapshot_start=argv[2];
  if(argc>3)
	snapshot_end=argv[3];
  else
	snapshot_end=snapshot_start;
}
using namespace std;
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "simulation_io.h"

void Snapshot_t::FormatSnapshotNumber(stringstream &ss)
{
  ss << setw(3) << setfill('0') << SnapshotIndex;
}

void Snapshot_t::Load(int snapshot_index)
{
  
}


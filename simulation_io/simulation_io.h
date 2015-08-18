#ifndef SIMULATION_IO_H_INCLUDED
#define SIMULATION_IO_H_INCLUDED

#include <iostream>
#include <sstream>
#include "../datatypes.h"

class Snapshot_t
{
public:
  int SnapshotIndex;
  List_t <Particle_t> ParticleList;
  Snapshot_t(int snapshot_index=0):ParticleList(),SnapshotIndex(snapshot_index)
  {
  }
  void Load(int snapshot_index);
  void FormatSnapshotNumber(std::stringstream &ss);
  void Clear()
  {
	ParticleList.Clear();
	SnapshotIndex=0;
  }
};

#endif
/* Each track records the life of a subhalo from birth time to final snapshot of the simulation */
#ifndef TRACKS_H_INCLUDED
#define TRACKS_H_INCLUDED

#include <iostream>
#include <new>

#include "datatypes.h"
#include "simulation_io/snapshot.h"
#include "simulation_io/halo.h"

class Track
{
public:
  HBTInt TrackId;
  List_t <SubHalo_t> SubHaloList;
};

#endif
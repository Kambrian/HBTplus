#ifndef GADGET_GROUP_IO_INCLUDED
#define GADGET_GROUP_IO_INCLUDED

#include "../halo.h"

namespace GadgetGroup
{
struct GroupV4Header_t
{
  int Ngroups;
  int Nsubgroups;
  int Nids;
  int TotNgroups;
  int TotNsubgroups;
  int TotNids;
  int num_files;//long? no, but padding may exist.
  double time;
  double redshift;
  double HubbleParam;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  int flag_doubleprecision;//long? no, but padding may exist
};

extern HBTInt Load(int SnapshotId, vector <Halo_t> &Halos);
extern bool IsGadgetGroup(const string &GroupFileFormat);

}

#endif
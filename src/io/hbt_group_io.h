/*
 * Copyright (C) 2007-2026 Jiaxin Han and contributors
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef HBT_GROUP_IO_HEADER_INCLUDED
#define HBT_GROUP_IO_HEADER_INCLUDED

#include "../datatypes.h"
#include "../halo.h"
namespace HBTGroupIO
{
  extern void GetHDFFileName(string &filename, int SnapshotIndex);
  extern HBTInt LoadHDFGroups(int snapshot_index, int snapshot_id, vector <Halo_t> &Halos);
  extern void SaveHDFGroups(int snapshot_index, int snapshot_id, vector <vector<HBTInt> > &halos);
}
#endif
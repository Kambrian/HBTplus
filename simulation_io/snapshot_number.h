#ifndef SNAPSHOT_NUMBER_H_INCLUDED
#define SNAPSHOT_NUMBER_H_INCLUDED

#include <iostream>
#include <sstream>
#include <iomanip>
#include <assert.h>
// #include <cstdlib>
// #include <cstdio>
// #include <unordered_map>
#include "../datatypes.h"
#include "../config_parser.h"

class SnapshotNumber_t
{
protected:  
  int SnapshotIndex;
  int SnapshotId;
public:
  SnapshotNumber_t()
  {
	SnapshotIndex=SpecialConst::NullSnapshotId;
	SnapshotId=SpecialConst::NullSnapshotId;
  }
  SnapshotNumber_t(SnapshotNumber_t & sn):SnapshotId(sn.SnapshotId), SnapshotIndex(sn.SnapshotIndex)
  {
  }
  SnapshotNumber_t & operator=(SnapshotNumber_t &sn)
  {
	SnapshotIndex=sn.SnapshotIndex;
	SnapshotId=sn.SnapshotId;
  }
  void ResetSnapshotNumber()
  {//reset is not destructon! when destructor is called, the data content no matter matters.
	SnapshotIndex=SpecialConst::NullSnapshotId;
	SnapshotId=SpecialConst::NullSnapshotId;
  }
  void FormatSnapshotId(std::stringstream &ss);
  void SetSnapshotIndex(int snapshot_index);
  int GetSnapshotIndex() const;
  int GetSnapshotId() const;
};
inline int SnapshotNumber_t::GetSnapshotIndex() const
{
  return SnapshotIndex;
}
inline int SnapshotNumber_t::GetSnapshotId() const
{
  return SnapshotId;
}
inline void SnapshotNumber_t::FormatSnapshotId(stringstream& ss)
{
  ss << std::setw(3) << std::setfill('0') << SnapshotId;
}
inline void SnapshotNumber_t::SetSnapshotIndex(int snapshot_index)
{
  assert(snapshot_index>=HBTConfig.MinSnapshotIndex&&snapshot_index<=HBTConfig.MaxSnapshotIndex);
//   assert(SpecialConst::NullSnapshotId!=snapshot_index);
  SnapshotIndex=snapshot_index; 
  if(HBTConfig.SnapshotIdList.empty())
	SnapshotId=SnapshotIndex;
  else
	SnapshotId=HBTConfig.SnapshotIdList[SnapshotIndex];
}

#endif
#include "mymath.h"
#include "linkedlist.h"

class PositionSample_t: public PositionData_t
{
private:
  int ThreadId, NumThreads;
  PositionData_t *Data;
  HBTInt np;
public:
  void init(int ithread, int nthread,  PositionData_t *data)
  {
    Data=data;
    ThreadId=ithread;
    NumThreads=nthread;
    HBTInt n0=Data->size();
    np=n0/nthread+((n0%nthread)>ithread);
  }
  const HBTxyz & operator [](HBTInt i) const
  {
    return (*Data)[restore_id(i)];
  }
  size_t size() const
  {
    return np;
  }
  HBTInt restore_id(HBTInt i) const
  {
    return i*NumThreads+ThreadId;
  }
  void restore_id(vector <HBTInt> &ids) const
  {
    for(auto &&i: ids)
      i=restore_id(i);
  }
};

class LinkedlistPara_t
{
private:
  vector<Linkedlist_t> LLs;
  vector <PositionSample_t> Samples;
public:
  LinkedlistPara_t(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false); 
  void SearchSphere(HBTReal radius, const HBTxyz &searchcenter, vector <HBTInt> &found_ids, int nmax_guess=8, HBTReal rmin=-1.);
  void SearchSphereSerial(HBTReal radius, const HBTxyz &searchcenter, vector <HBTInt> &found_ids, int nmax_guess=8, HBTReal rmin=-1.);
};



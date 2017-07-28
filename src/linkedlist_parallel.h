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
  void restore_id(vector <LocatedParticle_t> &particles) const
  {
    for(auto &&p: particles)
      p.id=restore_id(p.id);
  }
};

class LinkedlistPara_t
{
private:
  vector<Linkedlist_t> LLs;
  vector <PositionSample_t> Samples;
public:
  LinkedlistPara_t(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false); 
  void SearchSphere(HBTReal radius, const HBTxyz &searchcenter, vector <LocatedParticle_t> &founds)
  {
    SearchShell(-1., radius, searchcenter, founds);
  }
  void SearchSphereSerial(HBTReal radius, const HBTxyz &searchcenter, vector <LocatedParticle_t> &founds)
  {
    SearchShellSerial(-1., radius, searchcenter, founds);
  }
  void SearchShell(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter, vector <LocatedParticle_t> &founds);
  void SearchShellSerial(HBTReal rmin, HBTReal rmax, const HBTxyz &searchcenter, vector <LocatedParticle_t> &founds);
};



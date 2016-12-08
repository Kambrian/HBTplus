#include "mymath.h"

//TODO:discard the fortran-style ll; use struct or indexed table to parallelize the linklist!
struct LocatedParticle_t
{
  HBTInt id;
  HBTReal d; //distance
  LocatedParticle_t()=default;
  LocatedParticle_t(HBTInt id, HBTReal d):id(id),d(d)
  {}
};
class PositionData_t
{
public:
  virtual const HBTxyz & operator [](HBTInt i) const=0;
  /*virtual const HBTReal GetPos(HBTInt i, int j) const
  {
    return (*this)[i][j];
  }*/
  virtual size_t size() const=0;
};
class Linkedlist_t
{
private:
  int NDiv, NDiv2;
  bool PeriodicBoundary;
  HBTReal BoxSize, BoxHalf;
  vector <HBTInt> HOC;
  vector <HBTInt> List;
  HBTReal Range[3][2];
  HBTReal Step[3];
  PositionData_t *Particles;
  int RoundGridId(int i);
  int ShiftGridId(int i);
  int FixGridId(int i);
  HBTInt Sub2Ind(int i, int j, int k);
  HBTInt GetHOC(int i, int j, int k);
  HBTInt GetHOCSafe(int i, int j, int k);
  HBTReal Distance(const HBTxyz &x, const HBTxyz &y);
public:
  Linkedlist_t()=default;
  Linkedlist_t(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false)
  {
    build(ndiv, data, boxsize, periodic);
  }
  void build(int ndiv, PositionData_t *data, HBTReal boxsize=0., bool periodic=false);
  void SearchSphere(HBTReal radius, const HBTxyz &searchcenter, vector <LocatedParticle_t> &founds, int nmax_guess=8, HBTReal rmin=-1.);
};

inline int Linkedlist_t::RoundGridId(int i)
//to correct for rounding error near boundary
{
  return i<0?0:(i>=NDiv?NDiv-1:i);
}
inline int Linkedlist_t::ShiftGridId(int i)
/*to correct for periodic conditions; 
only applicable when def PERIODIC_BDR and ll.UseFullBox=1 */
{
      i=i%NDiv;
      if(i<0) i+=NDiv;
      return i;
}
inline int Linkedlist_t::FixGridId(int i)
{
  if(PeriodicBoundary)
    return ShiftGridId(i);
  return RoundGridId(i);
}
inline HBTInt Linkedlist_t::Sub2Ind(int i, int j, int k)
{
  return i+j*NDiv+k*NDiv2;
}
inline HBTInt Linkedlist_t::GetHOC(int i, int j, int k)
{
  return HOC[Sub2Ind(i,j,k)];
}
inline HBTInt Linkedlist_t::GetHOCSafe(int i, int j, int k)
{
  return HOC[Sub2Ind(FixGridId(i), FixGridId(j), FixGridId(k))];
}
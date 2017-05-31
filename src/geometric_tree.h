#ifndef GEOMETRIC_TREE_H_INCLUDED
#define GEOMETRIC_TREE_H_INCLUDED
/*specialized octtree for spatial searching*/

#include "oct_tree.h"

typedef TreeCell_t<HBTInt> GeoTreeCell_t;

class GeoTree_t: public OctTree_t<GeoTreeCell_t>
{
private:
  void ProcessNode(HBTInt nodeid, HBTInt nextid, int sonid, HBTInt &mass, double len, const double center[3]);
  void FillNodeCenter(HBTInt nodeid, const double center[3]);
  void UpdateInternalNodes(HBTInt no,HBTInt sib,double len, const double center[3]);
public:
  void Search(const HBTxyz &searchcenter, HBTReal radius, vector <LocatedParticle_t> &founds);
  HBTInt NearestNeighbour(const HBTxyz &searchcenter, HBTReal rguess);
  double SphDensity(const HBTxyz &cen, HBTReal & rguess);
};

#endif	















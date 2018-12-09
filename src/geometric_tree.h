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
  int NumNeighbourSPH;
public:
  GeoTree_t():OctTree_t<GeoTreeCell_t>(), NumNeighbourSPH(64)
  {
  }
  void Search(const HBTxyz &searchcenter, HBTReal radius, ParticleCollector_t &collector);
  HBTInt NearestNeighbour(const HBTxyz &searchcenter, HBTReal rguess);
  double SphDensity(const HBTxyz &cen, HBTReal & rguess);
  int GetNumNeighbourSPH()
  {
    return NumNeighbourSPH;
  }
  int SetNumNeighbourSPH(int num_neighbour)
  {
    NumNeighbourSPH=num_neighbour;
  }
};

inline HBTReal GuessNeighbourRange(HBTInt n_neighbours, HBTReal number_density_guess)
{
  return pow(3 * n_neighbours / (4 * 3.141593) / number_density_guess, 1.0 / 3);
}

#endif	















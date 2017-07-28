#ifndef GRAVITY_TREE_H_INCLUDED
#define GRAVITY_TREE_H_INCLUDED
/*specialized octtree for gravity calculation*/

#include "oct_tree.h"

typedef TreeCell_t<HBTReal> GravityTreeCell_t;
class GravityTree_t:public OctTree_t<GravityTreeCell_t>
{
private:
  void ProcessNode(HBTInt nodeid, HBTInt nextid, int sonid, double &mass, double CoM[3], double len, const double center[3]); 
  void FillNodeCenter(HBTInt nodeid, const double center[3], double CoM[3], double mass);
  void UpdateInternalNodes(HBTInt no,HBTInt sib,double len, const double center[3]);
public:
  double EvaluatePotential(const HBTxyz &targetPos, const HBTReal targetMass=0.);
  double BindingEnergy(const HBTxyz &targetPos, const HBTxyz &targetVel, const HBTxyz &refPos, const HBTxyz &refVel, const HBTReal targetMass=0.);
};

#endif	















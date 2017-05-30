//ToDo: separate GravityTree and SpatialTree class; add searching and density functions; neareast neighbour search by inserting into the nodes?
#ifndef FOF_BUILDER_H_INCLUDED
#define FOF_BUILDER_H_INCLUDED
#include <exception>

#include "gravity_tree.h"

class FoFBuilder_t: public OctTree_t
{
private:
  double LinkLength2, LinkLengthNode;
  bool IsPeriodic;
public:
  vector <HBTInt> GrpLen, GrpTags;
  HBTReal LinkLength;
  FoFBuilder_t(HBTReal linklength, const Snapshot_t &snapshot): LinkLength(linklength), OctTree_t()
  {
    IsPeriodic=HBTConfig.PeriodicBoundaryOn;
    
    LinkLength2=LinkLength*LinkLength;
    LinkLengthNode=LinkLength/3.;//check for node entirely when nodesize is small enough
    cout<<"Building tree...\n"<<flush;
    Build(snapshot, 0, false);
    GrpTags.assign(snapshot.size(), -1);
  }
  void Link();
  HBTInt TagFriendsOfFriends(HBTInt seed, HBTInt grpid);
  void PullFriends(OctTree_t::OctTreeCell_t& rootnode, const HBTxyz& searchcenter, HBTInt grpid, vector< HBTInt >& friends);
  void PullNode(OctTreeCell_t &node, HBTInt grpid, vector <HBTInt> &friends);
};

#endif	















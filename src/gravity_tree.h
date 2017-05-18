/* this oct-tree code is adopted from SUBFIND and modified for HBT */
//ToDo: separate GravityTree and SpatialTree class; add searching and density functions; neareast neighbour search by inserting into the nodes?
#ifndef TREE_H_INCLUDED
#define TREE_H_INCLUDED
#include <exception>

#include "datatypes.h"
#include "snapshot.h"

class OctTreeExceeded_t : public exception
{
private:
  string msg;
public:
  OctTreeExceeded_t(const string & message)
  {
	msg=message;
  }
  const char * what () const throw ()
  {
	return msg.c_str();
  }
  ~OctTreeExceeded_t() throw()
  {}
};
class OctTree_t
{
private:  
  union OctTreeCell_t
  {
	HBTInt sons[8];		/*!< temporary pointers to daughter nodes */
	struct
	{
	HBTReal s[3];               /*!< center of mass of node (gravity tree); geocenter for geotree*/
	HBTReal len;		/*!< sidelength of treenode */
	HBTReal mass;            /*!< mass of node (gravity tree); counts of particles for geotree */
	HBTInt sibling;         /*!< this gives the next node in the walk in case the current node can be used */
	HBTInt nextnode;        /*!< this gives the next node in case the current node needs to be opened */
	}way;
  };
  /*the storage*/
  OctTreeCell_t *Cells; 
  HBTInt * NextnodeFromParticle; /* next node for each particle. Particles are the first NumPart nodes, and cells are the remaining nodes.*/
  /*auxilliary pointer*/
  OctTreeCell_t *Nodes;   /* =Cells-NumberOfParticles. the nodes are labelled from 0 to NumPart+NumNodes-1, so that nodeid=0~NumPart-1 are particles, and nodeid>=NumPart are cells */
  size_t MaxNumberOfCells, MaxNumberOfParticles;
  HBTInt MaxNodeId;
  const Snapshot_t * Snapshot;
  HBTInt NumberOfParticles; //alias to Snapshot->GetSize().
  HBTInt & RootNodeId; //alias to NumberOfParticles
  void UpdateInternalNodes(HBTInt no,HBTInt sib,double len, const double center[3]);
  void UpdateSubnode(HBTInt son, HBTInt sib, double len, const double center[3], int subindex);
public:
  bool IsGravityTree;
  OctTree_t(): MaxNumberOfCells(0), MaxNumberOfParticles(0), MaxNodeId(0), NumberOfParticles(0), RootNodeId(NumberOfParticles), IsGravityTree(true)
  {
  }
  void Reserve(const size_t max_num_part);
  HBTInt Build(const Snapshot_t &snapshot, HBTInt num_part=0, bool ForGravity=true);
  void Clear();
  double EvaluatePotential(const HBTxyz &targetPos, const HBTReal targetMass=0.);
  double BindingEnergy(const HBTxyz &targetPos, const HBTxyz &targetVel, const HBTxyz &refPos, const HBTxyz &refVel, const HBTReal targetMass=0.);
  void Search(const HBTxyz &searchcenter, HBTReal radius, vector <LocatedParticle_t> &founds);
  HBTInt NearestNeighbour(const HBTxyz &searchcenter, HBTReal rguess);
  double SphDensity(const HBTxyz &cen, HBTReal & rguess);
  HBTInt TagFriendsOfFriends(HBTInt seed, HBTInt grpid, vector <HBTInt> &GroupTags, double LinkLength);
  void TagNode(OctTreeCell_t & node, HBTInt grpid, vector <HBTInt> &GroupTags, vector <HBTInt> &friends);
  ~OctTree_t()
  {
	Clear();
  }
};

/*
class GravityTree_t: public OctTree_t
{
  //you can override the UpdateInternalNodes function (set to virtual in base class first).
}
*/

extern void treesearch_linkgrp(HBTReal radius, const Snapshot_t &snapshot, vector <HBTInt> &GrpLen, vector <HBTInt> &GrpTags);

#endif	















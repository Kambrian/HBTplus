/* this oct-tree code is adopted from SUBFIND with minor modifications for HBT */
#ifndef TREE_H_INCLUDED
#define TREE_H_INCLUDED
#include <exception>

#include "../datatypes.h"
#include "../simulation_io/snapshot.h"

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
	HBTReal s[3];               /*!< center of mass of node */
	HBTReal len;		/*!< sidelength of treenode */
	HBTReal mass;            /*!< mass of node */
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
  HBTInt NumberOfParticles, MaxNodeId;
  const HBTInt *ParticleList;
  const Snapshot_t * Snapshot;
  void UpdateInternalNodes(HBTInt no,HBTInt sib,double len);
public:
  OctTree_t(): MaxNumberOfCells(0), MaxNumberOfParticles(0), MaxNodeId(0), NumberOfParticles(0)
  {
  }
  void Reserve(const size_t max_num_part);
  HBTInt Build(const HBTInt num_part, const HBTInt * particles, const Snapshot_t &snapshot);
  void Clear();
  double EvaluatePotential(const HBTReal targetPos[3], const HBTReal targetMass=0.);
  double BindingEnergy(const HBTxyz &targetPos, const HBTxyz &targetVel, const HBTxyz &refPos, const HBTxyz &refVel, const HBTReal targetMass=0.);
  ~OctTree_t()
  {
	Clear();
  }
};


#endif	















/* this oct-tree code is adopted from SUBFIND with minor modifications for HBT */
#ifndef TREE_H_INCLUDED
#define TREE_H_INCLUDED

#include "../datatypes.h"

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
	HBTInt mass;            /*!< mass of node */
	HBTInt sibling;         /*!< this gives the next node in the walk in case the current node can be used */
	HBTInt nextnode;        /*!< this gives the next node in case the current node needs to be opened */
	}way;
  };
  /*the storage*/
  OctTreeCell_t *Cells; 
  HBTInt * NextnodeFromParticle; /* next node for each particle. Particles are the first NumPart nodes, and cells are the remaining nodes.*/
  /*auxilliary pointer*/
  OctTreeCell_t *Nodes;   /*!< this is a shited pointer to access the cells with nodeid, such that Nodes[NumPart] gives the first cell; the nodes are labelled from 0 to NumPart+NumNodes-1, so that nodeid=0~NumPart-1 are particles, and nodeid>=NumPart are cells */
  size_t MaxNumberOfCells, MaxNumberOfNodes;
  HBTInt NumberOfParticles;
  HBTInt *ParticleList;
  HBTxyz *ParticlePosition;
  void UpdateInternalNodes(HBTInt no,HBTInt sib,double len);
  void Allocate(size_t num_part);
public:
  OctTree_t(): MaxNumberOfCells(0), MaxNumberOfNodes(0), NumberOfParticles(0)
  {
  }
  HBTInt Build(HBTInt num_part, HBTInt * particles, HBTReal part_pos[][3]);
  void Clear();
  double EvaluatePotential(HBTReal targetPos[3]);
  ~OctTree_t()
  {
	Clear();
  }
};

#endif	















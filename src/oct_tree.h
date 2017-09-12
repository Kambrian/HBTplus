#ifndef TREE_H_INCLUDED
#define TREE_H_INCLUDED
/* template for OctTree; specialize into GravityTree_t and GeometricTree_t, depending on whether the center of mass or geometric center is recorded.*/

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

template <class T> union TreeCell_t
{
  typedef T MassType_t;
  HBTInt sons[8];		/*!< temporary pointers to daughter nodes */
  struct
  {
    HBTReal s[3];               /*!< center of mass of node (gravity tree); geocenter for geotree*/
    HBTReal len;		/*!< sidelength of treenode */
    T mass;            /*!< mass of node (gravity tree); counts of particles for geotree */
    HBTInt sibling;         /*!< this gives the next node in the walk in case the current node can be used */
    HBTInt nextnode;        /*!< this gives the next node in case the current node needs to be opened */
  }way;
  TreeCell_t(){};
  TreeCell_t(HBTInt i): sons{i}
  {
  }
};

template <class CellT> 
class OctTree_t
{
protected:  
  typedef CellT OctTreeCell_t;
  /*the storage*/
  vector <OctTreeCell_t> Cells;
  OctTreeCell_t *Nodes;   /* =Cells-NumberOfParticles. the nodes are labelled from 0 to NumPart+NumNodes-1, so that nodeid=0~NumPart-1 are particles, and nodeid>=NumPart are cells */
  vector <HBTInt> NextnodeFromParticle; /* next node for each particle. Particles are the first NumPart nodes, and cells are the remaining nodes.*/
  const Snapshot_t * Snapshot;
  HBTInt NumberOfParticles; //alias to Snapshot->GetSize().
  HBTInt & RootNodeId; //alias to NumberOfParticles
private:  
  virtual void UpdateInternalNodes(HBTInt no,HBTInt sib,double len, const double center[3])=0;
public:
  OctTree_t(): NumberOfParticles(0), RootNodeId(NumberOfParticles)
  {
  }
  void Reserve(const size_t max_num_part);
  HBTInt Build(const Snapshot_t &snapshot, HBTInt num_part=0);
  void AppendCell();
  void Clear();
  ~OctTree_t()
  {
	Clear();
  }
};

#include "oct_tree.tpp"

#endif	















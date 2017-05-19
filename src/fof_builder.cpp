#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "fof_builder.h"

void FoFBuilder_t::Link()
/* link particles in the given snapshot into groups.
 * Output: filled GrpLen and GrpTags (0~Ngroups-1), down to mass=1 (diffuse particles)
 * */
{
cout<<"Linking Groups...\n"<<flush;
HBTInt grpid=0;
HBTInt printstep=Snapshot->size()/100, progress=printstep;
cout<<"00%"<<flush;
for(HBTInt i=0;i<Snapshot->size();i++)
{
	if(GrpTags[i]<0)
	{
		if(i>=progress)
		{
		  cout<<"\b\b\b"<<setw(2)<<progress/printstep<<"%"<<flush;
		  progress+=printstep;
		}
		GrpTags[i]=grpid; //infect the seed particle
		HBTInt grplen=1+TagFriendsOfFriends(i,grpid);
		GrpLen.push_back(grplen);
		grpid++;
	}
}

cout<<"Found "<<GrpLen.size()<<" Groups\n";
}

inline void FoFBuilder_t::PullNode(OctTreeCell_t &node, HBTInt grpid, vector <HBTInt> &friends)
/*tag entire node with grpid and clear it*/
{
  HBTInt end_node_id=node.way.sibling;
  HBTInt node_id=node.way.nextnode;
  while(node_id!=end_node_id)
  {
    if(node_id<NumberOfParticles)
    {
      if(GrpTags[node_id]<0)
      {
	GrpTags[node_id]=grpid;
	friends.push_back(node_id);
      }
      node_id=NextnodeFromParticle[node_id];
    }
    else
      node_id=Nodes[node_id].way.nextnode;//open it straight-away
  }
  //short-circuit this node
  node.way.nextnode=node.way.sibling;
  node.way.mass=0.;
}

void FoFBuilder_t::PullFriends(OctTreeCell_t &rootnode, const HBTxyz &searchcenter, HBTInt grpid, vector <HBTInt> &friends)
/*tag all the particles inside node that are directly linked to seed with grpid.
   * this function should only be used with a geometric tree (ie, IsGravityTree=false)
   **/
{
  double x0=searchcenter[0], y0=searchcenter[1], z0=searchcenter[2];
  HBTInt node_id = rootnode.way.nextnode, end_node_id=rootnode.way.sibling;
  
  HBTInt nold=friends.size();
  while(node_id !=end_node_id)
  {
    if(node_id < NumberOfParticles)		/* single particle */
    {
      HBTInt pid=node_id;
      node_id = NextnodeFromParticle[pid];
      if(GrpTags[pid]>=0) //already tagged 
	continue;  
      
      auto &pos=Snapshot->GetComovingPosition(pid);
      double dx = pos[0] - x0;
      if(IsPeriodic) dx=NEAREST(dx);
      if(dx >LinkLength || dx <-LinkLength)	continue;
      
      double dy = pos[1] - y0;
      if(IsPeriodic) dy=NEAREST(dy);
      if(dy >LinkLength || dy <-LinkLength)	continue;
      
      double dz = pos[2] - z0;
      if(IsPeriodic) dz=NEAREST(dz);
      if(dz >LinkLength || dz <-LinkLength)	continue;
      
      double r2 = dx * dx + dy * dy + dz * dz;
      
      if(r2 < LinkLength2) //confirm the infection (fill the stack)
      {
	GrpTags[pid]=grpid;
	friends.push_back(pid);
      }
    }
    else
    {
      auto & node = Nodes[node_id];
      node_id = node.way.sibling;	/* in case the node can be discarded */
      if(node_id==node.way.nextnode)//closed node
	continue;
      
      double lenhalf=node.way.len/2.;
      double rmax=lenhalf+LinkLength;
      
      auto &pos=node.way.s;
      double dx = pos[0] - x0;
      if(IsPeriodic) dx=NEAREST(dx);
      if(dx > rmax || dx <-rmax)	continue;
      
      double dy = pos[1] - y0;
      if(IsPeriodic) dy=NEAREST(dy);
      if(dy > rmax || dy <-rmax)	continue;
      
      double dz = pos[2] - z0;
      if(IsPeriodic) dz=NEAREST(dz);
      if(dz > rmax || dz <-rmax)	continue;
      
      /*if(node.way.mass>3&&lenhalf<LinkLengthNode)//a small node (max lenhalf=LinkLength/sqrt(3)), check if it fits entirely
      {
	dx=dx>0?lenhalf+dx:lenhalf-dx;
	dy=dy>0?lenhalf+dy:lenhalf-dy;
	dz=dz>0?lenhalf+dz:lenhalf-dz;
	double r2 = dx*dx+dy*dy+dz*dz;
	if(r2 < LinkLength2) //entire node can be used
	{
	  PullNode(node, grpid, friends);
	  continue;
	}
      }*/
      
      PullFriends(node, searchcenter, grpid, friends);
    }
  }
  
  rootnode.way.mass-=(friends.size()-nold);
  if(rootnode.way.mass==0.)
    rootnode.way.nextnode=end_node_id;
}

HBTInt FoFBuilder_t::TagFriendsOfFriends(HBTInt seed, HBTInt grpid)
/*tag all the particles that are linked (both directly and indirectly) to seed with grpid
   * Note if system stack size is too small, this recursive routine may crash
   * in that case you should set:  ulimit -s unlimited  (bash) before running.
   * this function should only be used with a geometric tree (ie, IsGravityTree=false)
   **/
{
  auto &searchcenter=Snapshot->GetComovingPosition(seed);
  
  vector <HBTInt> friends;
  PullFriends(Nodes[RootNodeId], searchcenter,  grpid, friends);
  
  HBTInt nfriends=friends.size();
  for(auto &&pid: friends)
    nfriends+=TagFriendsOfFriends(pid, grpid);
  
  return nfriends; //excluding the seed particle
}

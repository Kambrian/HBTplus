#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "datatypes.h"
#include "gravity_tree.h"

HBTReal GuessNeighbourRange(HBTInt n_neighbours, HBTReal number_density_guess)
{
  return pow(3 * n_neighbours / (4 * 3.141593) / number_density_guess, 1.0 / 3);
}

HBTInt OctTree_t::NearestNeighbour(const HBTxyz & cen, HBTReal rguess)
{
  vector <LocatedParticle_t> founds;
  Search(cen, rguess, founds);
  while(founds.empty())         //WARNING: dead loop if tree is empty.
  {
    rguess *= 1.26;//double the guess volume
    Search(cen, rguess, founds);
  }
  return min_element(founds.begin(), founds.end(), CompLocatedDistance)->id;	
}

void OctTree_t::Search(const HBTxyz & searchcenter, HBTReal radius, vector <LocatedParticle_t> &founds)
{/*find a list of particles from the tree, located within radius around searchcenter,
  * and APPEND their particle_id and distance^2 to founds */
  bool IsPeriodic=HBTConfig.PeriodicBoundaryOn;
  double x0=searchcenter[0], y0=searchcenter[1], z0=searchcenter[2];
  double h2 = radius * radius;

  HBTInt numngb = 0;
  HBTInt node_id = RootNodeId;

  while(node_id >= 0)
    {
      if(node_id < RootNodeId)		/* single particle */
	{
	  HBTInt pid=node_id;	  
	  node_id = NextnodeFromParticle[node_id];
	  
	  auto &pos=Snapshot->GetComovingPosition(pid);
	  double dx = pos[0] - x0;
	  if(IsPeriodic) dx=NEAREST(dx);
	  if(dx < -radius || dx > radius)
	    continue;
	  
	  double dy = pos[1] - y0;
	  if(IsPeriodic) dy=NEAREST(dy);
	  if(dy < -radius || dy > radius)
	    continue;

	  double dz = pos[2] - z0;
	  if(IsPeriodic) dz=NEAREST(dz);
	  if(dz < -radius || dz > radius)
	    continue;

	  double r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	      founds.emplace_back(Snapshot->GetMemberId(pid), r2);
	}
      else
	{
	  auto &node = Nodes[node_id];

	  node_id = node.way.sibling;	/* in case the node can be discarded */
	  double rmax=node.way.len;
	  if(!IsGravityTree) rmax/=2;
	  rmax+=radius;
	  
	  auto &pos=node.way.s;
	  double dx = pos[0] - x0;
	  if(IsPeriodic) dx=NEAREST(dx);
	  if(dx < -rmax || dx > rmax)
	    continue;
	  
	  double dy = pos[1] - y0;
	  if(IsPeriodic) dy=NEAREST(dy);
	  if(dy < -rmax || dy > rmax)
	    continue;

	  double dz = pos[2] - z0;
	  if(IsPeriodic) dz=NEAREST(dz);
	  if(dz < -rmax || dz > rmax)
	    continue;

	  node_id = node.way.nextnode;	/* ok, we need to open the node */
	}
    }
}

#define SPH_DENS_NGB 64
double OctTree_t::SphDensity(const HBTxyz &cen, HBTReal & hguess)
{
  vector <LocatedParticle_t> founds;
  Search(cen, hguess, founds);
  int numngb=founds.size();
  while(numngb<SPH_DENS_NGB)
  {
    if(numngb)	
      hguess *= pow(1.*SPH_DENS_NGB/numngb,1.0/3.0)*1.1;//update hguess adaptively, and conservatively to keep it slightly larger
    else  //zero ngb, double hguess
	hguess *= 2.;
      
      founds.clear();
    Search(cen, hguess, founds);
    numngb=founds.size();
  }
  
  auto pivot_particle=founds.begin()+SPH_DENS_NGB-1;
  nth_element(founds.begin(), founds.end(), pivot_particle, CompLocatedDistance);
  double h=sqrt(pivot_particle->d);
  // 	h=sqrtf(h);
  hguess=h*1.01;
  double hinv3 = 1.0 / (h * h * h);
  
  double rho=0.;
  for(auto it=founds.begin(); it <=pivot_particle; ++it)
  {
    double r = sqrt(it->d);
    double u = r / h, wk;
    
    if(u < 0.5)
      wk = hinv3 * (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
    else
      wk = hinv3 * 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
    
    rho += wk;
  }
  return rho;
}

static double LinkLength, LinkLength2;
#pragma omp threadprivate(LinkLength,LinkLength2)

void treesearch_linkgrp(HBTReal radius, const Snapshot_t &snapshot, vector <HBTInt> &GrpLen, vector <HBTInt> &GrpTags)
/* link particles in the given snapshot into groups.
 * Output: filled GrpLen and GrpTags (0~Ngroups-1), down to mass=1 (diffuse particles)
 * */
{
GrpTags.assign(snapshot.size(), -1);
LinkLength=radius;
LinkLength2=radius*radius;

cout<<"Building tree...\n"<<flush;
OctTree_t tree;
tree.Build(snapshot, 0, false);

cout<<"Linking Groups...\n"<<flush;
HBTInt grpid=0;
HBTInt printstep=snapshot.size()/100, progress=printstep;
cout<<"00%%"<<flush;
for(HBTInt i=0;i<snapshot.size();i++)
{
	if(GrpTags[i]<0)
	{
		if(i>=progress)
		{
		  cout<<"\b\b\b"<<setw(2)<<progress/printstep<<"%%"<<flush;
		  progress+=printstep;
		}
		GrpTags[i]=grpid; //infect the seed particle
		HBTInt grplen=1+tree.TagFriendsOfFriends(i,grpid, GrpTags);
		GrpLen.push_back(grplen);
		grpid++;
	}
}

cout<<"Found "<<GrpLen.size()<<" Groups\n";
}

HBTInt OctTree_t::TagFriendsOfFriends(HBTInt seed, HBTInt grpid,
				      vector <HBTInt> &group_tags)
/*tag all the particles that are linked to seed with grpid
   * Note if system stack size is too small, this recursive routine may crash
   * in that case you should set:  ulimit -s unlimited  (bash) before running.
   **/
{
  bool IsPeriodic=HBTConfig.PeriodicBoundaryOn;
  auto &searchcenter=Snapshot->GetComovingPosition(seed);
  double x0=searchcenter[0], y0=searchcenter[1], z0=searchcenter[2];
  
  HBTInt node_id = RootNodeId;
  
  vector <HBTInt> friends;
  while(node_id >= 0)//infect neighbours
  {
    if(node_id < NumberOfParticles)		/* single particle */
    {
      HBTInt pid=node_id;
      node_id = NextnodeFromParticle[pid];
      if(group_tags[pid]>=0) //already tagged 
	continue;  
      
      auto &pos=Snapshot->GetComovingPosition(pid);
      double dx = pos[0] - x0;
      if(IsPeriodic) dx=NEAREST(dx);
      if(dx < -LinkLength || dx > LinkLength)
	continue;
      
      double dy = pos[1] - y0;
      if(IsPeriodic) dy=NEAREST(dy);
      if(dy < -LinkLength || dy > LinkLength)
	continue;
      
      double dz = pos[2] - z0;
      if(IsPeriodic) dz=NEAREST(dz);
      if(dz < -LinkLength || dz > LinkLength)
	continue;
      
      double r2 = dx * dx + dy * dy + dz * dz;
      
      if(r2 < LinkLength2) //confirm the infection (fill the stack)
      {
	group_tags[pid]=grpid;
	friends.push_back(pid);
      }
    }
    else
    {
      auto & node = Nodes[node_id];
      
      node_id = node.way.sibling;	/* in case the node can be discarded */
      double rmax=node.way.len;
      if(!IsGravityTree) rmax/=2;
      rmax+=LinkLength;
      
      auto &pos=node.way.s;
      double dx = pos[0] - x0;
      if(IsPeriodic) dx=NEAREST(dx);
      if(dx < -rmax || dx > rmax)
	continue;
      
      double dy = pos[1] - y0;
      if(IsPeriodic) dy=NEAREST(dy);
      if(dy < -rmax || dy > rmax)
	continue;
      
      double dz = pos[2] - z0;
      if(IsPeriodic) dz=NEAREST(dz);
      if(dz < -rmax || dz > rmax)
	continue;
      
      node_id = node.way.nextnode;	/* ok, we need to open the node */
    }
  }
  //TODO: optimize by deleting seed from tree? current tree structure not suitable for this
  
  HBTInt nfriends=friends.size();
  for(auto &&pid: friends)
    nfriends+=TagFriendsOfFriends(pid, grpid, group_tags);
  
  return nfriends; //excluding the seed particle
}

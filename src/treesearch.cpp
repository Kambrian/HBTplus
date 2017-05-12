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

HBTInt OctTree_t::NearestNeighbour(const HBTxyz & cen[3],HBTReal rguess)
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

void OctTree_t::Search(const HBTxyz & searchcenter, HBTReal radius, vector <LocateParticle_t> &founds)
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
	      founds.emplace_back(Snapshot->GetMemberId(pid), sqrt(r2));
	}
      else
	{
	  auto &node = &Nodes[node_id];

	  node_id = Nodes[node_id].way.sibling;	/* in case the node can be discarded */
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
  HBTInt i, n;
  double h, hinv3, wk, u, r, rho;
  
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
  h=pivot_particle->r;
  // 	h=sqrtf(h);
  hguess=h*1.01;
  hinv3 = 1.0 / (h * h * h);
  
  for(auto it=founds.begin(), rho = 0; it <=pivot_particle; ++it)
  {
    r = it->d;
    u = r / h;
    
    if(u < 0.5)
      wk = hinv3 * (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
    else
      wk = hinv3 * 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
    
    rho += wk;
  }
  return rho;
}

static HBTInt *InfectionStack, *pInfected;
static double LinkLength, LinkLength2;
#pragma omp threadprivate(InfectionStack,pInfected,LinkLength,LinkLength2)

HBTInt treesearch_linkgrp(HBTReal radius, HBTInt PIndex[], struct GroupData *GrpData)
/* To link DM particles into groups
 * Input: radius: link radius
 * 		  PIndex[]: list of particles to link
 *        grpdata: pointer to GroupData structure.
 * Output: filled grpdata.
 * Return value: number of groups found, down to mass=1 (diffuse particles)
 * */
{
HBTInt i,grpid;

//==init
GrpData->group_tags=mymalloc(sizeof(struct ParticleGroup)*GrpData->Np);
for(i=0;i<GrpData->Np;i++) 
{
	GrpData->group_tags[i].PIndex=PIndex[i];
	GrpData->group_tags[i].GrpID=-1;
}
GrpData->GrpLen=mymalloc(sizeof(HBTInt)*GrpData->Np);	

InfectionStack=mymalloc(sizeof(HBTInt)*GrpData->Np);//at most Np particles would be infected
pInfected=InfectionStack;//initial position
LinkLength=radius;
LinkLength2=radius*radius;

fprintf(logfile,"constructing tree...\n");fflush(logfile);
tree_tree_allocate(TREE_ALLOC_FACTOR*(size_t)GrpData->Np,GrpData->Np);
maketree(GrpData->Np,PIndex,Pdat.Pos);
	
fprintf(logfile,"Linking Groups...\n");fflush(logfile);
//~ printf("P%08d PhysicalConst::G%08d",0,0);
grpid=0;
HBTInt j,jj=1; j=NumPart/100;
fprintf(logfile,"00%%");fflush(logfile);
for(i=0;i<NumPart;i++)
{
	if(GrpData->group_tags[i].GrpID<0)
	{
		//~ int j;
		//~ for(j=0;j<19;j++) printf("\b");
		//~ printf("P%08d PhysicalConst::G%08d",i,grpid);fflush(stdout);
		if(i>=j){fprintf(logfile,"\b\b\b%02d%%",(int)jj);fflush(logfile);jj++;j=(NumPart/100.)*jj;}
		GrpData->group_tags[i].GrpID=grpid; //infect the seed particle
		GrpData->GrpLen[grpid]=1+treesearch_infect_particles(i,grpid, GrpData->group_tags, Pdat.Pos);
		grpid++;
		
	}
}

tree_tree_free();

GrpData->Ngrp=grpid;
GrpData->GrpLen=realloc(GrpData->GrpLen,sizeof(HBTInt)*grpid);
fprintf(logfile,"Found "HBTIFMT" Groups\n",grpid);fflush(logfile);

if(pInfected!=InfectionStack)
{
	fprintf(logfile,"Error: InfectionStack not cleaned in treesearch_linkgrp, current pos=%d\n",
						(int)(pInfected-InfectionStack));
}
myfree(InfectionStack);

return grpid;	
}

HBTInt OctTree_t::InfectParticles(HBTInt seed, HBTInt grpid,
		ParticleGroup_t &group_tags)
{
/*tag all the particles that are linked to seed with grpid
 * Note if system stack size is too small, this recursive routine may crash
 * in that case you should set:  ulimit -s unlimited  (bash) before running.
**/
 HBTInt numngb, totnumngb, p;
  double dx, dy, dz, r2, h2;
  union NODE *this;
  HBTReal *searchcenter;

  searchcenter=Snapshot->GetComovingPosition(group_tags[seed].ParticleId);
  //~ h2 = radius * radius;

  numngb = 0; 
  HBTInt node_id = RootNodeId;

  while(node_id >= 0)//infect neighbours
    {
      if(node_id < NumberOfParticles)		/* single particle */
	  {
		  *pInfected = node_id;
		  node_id = NextnodeFromParticle[*pInfected];
		  if(group_tags[*pInfected].GrpID>=0) //already tagged
			continue;  
		  
		  p = group_tags[*pInfected].PIndex;

		  dx = PPos[p][0] - searchcenter[0];
		  #ifdef PERIODIC_BDR
		  dx=NEAREST(dx);
		  #endif
		  if(dx < -LinkLength)
			continue;
		  if(dx > LinkLength)
			continue;

		  dy =PPos[p][1] - searchcenter[1];
		  #ifdef PERIODIC_BDR
		  dy=NEAREST(dy);
		  #endif
		  if(dy < -LinkLength)
			continue;
		  if(dy > LinkLength)
			continue;

		  dz = PPos[p][2] - searchcenter[2];
		  #ifdef PERIODIC_BDR
		  dz=NEAREST(dz);
		  #endif
		  if(dz < -LinkLength)
			continue;
		  if(dz > LinkLength)
			continue;

		  r2 = dx * dx + dy * dy + dz * dz;

		  if(r2 < LinkLength2) //confirm the infection (fill the stack)
			{
			  group_tags[*pInfected].GrpID=grpid;
			  numngb++;
			  pInfected++;
			}
	  }
      else
	  {
	  this = &Nodes[node_id];

	  node_id = Nodes[node_id].way.sibling;	/* in case the node can be discarded */
	//since way.s[3] is CoM rather than center of cube,compare with Len rather than Lenhalf to allow for misaligned CoM
	  if((NEAREST(this->way.s[0] - searchcenter[0]) + this->way.len) < -LinkLength)
	    continue;
	  if((NEAREST(this->way.s[0] - searchcenter[0]) - this->way.len) > LinkLength)
	    continue;
	  if((NEAREST(this->way.s[1] - searchcenter[1]) + this->way.len) < -LinkLength)
	    continue;
	  if((NEAREST(this->way.s[1] - searchcenter[1]) - this->way.len) > LinkLength)
	    continue;
	  if((NEAREST(this->way.s[2] - searchcenter[2]) + this->way.len) < -LinkLength)
	    continue;
	  if((NEAREST(this->way.s[2] - searchcenter[2]) - this->way.len) > LinkLength)
	    continue;

	  node_id = this->way.nextnode;	/* ok, we need to open the node */
	  }
    }
	
	//~ printf("ngb=%d ",numngb);
	totnumngb=numngb;
	while(numngb>0)//pop the stack
	{
		pInfected--;
		totnumngb+=treesearch_infect_particles(*pInfected,grpid, group_tags, PPos);
		numngb--;
	}
	return totnumngb; //total number of infected particles (excluding the seed particle)
}

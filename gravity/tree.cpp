/* this oct-tree code is adopted from SUBFIND with minor modifications for HBT */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "../mymath.h"
#include "tree.h"
#include "../config_parser.h"

void OctTree_t::UpdateInternalNodes(HBTInt no, HBTInt sib, double len)
{
	HBTInt j,jj,p,pp,mass,sons[8];
	double s[3];
	
		mass=0;
		s[0]=0;
		s[1]=0;
		s[2]=0;
		for(j=0;j<8;j++)
			sons[j]=Nodes[no].sons[j];//backup sons
		Nodes[no].way.len=len;
		Nodes[no].way.sibling=sib;
		for(j=0;sons[j]<0;j++);//find first son
		pp=sons[j];
		Nodes[no].way.nextnode=pp;
		for(jj=j+1;jj<8;jj++)//find sons in pairs,ie. find sibling
		{
			if(sons[jj]>=0)//ok, found a sibling
			{
				p=pp;
				pp=sons[jj];
				if(p<NumberOfParticles)
				{
					mass++;
					s[0]+=ParticlePosition[ParticleList[p]][0];
					s[1]+=ParticlePosition[ParticleList[p]][1];
					s[2]+=ParticlePosition[ParticleList[p]][2];
					NextnodeFromParticle[p]=pp;
				}
				else
				{
				if(len>=HBTConfig.TreeNodeResolution) 
					UpdateInternalNodes(p,pp,0.5*len);//only divide if above resolution; otherwise 
				//we didn't divide the node seriouly so we don't have finer node length
				else
					UpdateInternalNodes(p,pp,len);//get internal node info
				mass+=Nodes[p].way.mass;
				s[0]+=Nodes[p].way.s[0]*Nodes[p].way.mass;
				s[1]+=Nodes[p].way.s[1]*Nodes[p].way.mass;
				s[2]+=Nodes[p].way.s[2]*Nodes[p].way.mass;
				}
			}
		}
		if(pp<NumberOfParticles)//the last son
		{			
			mass++;
			s[0]+=ParticlePosition[ParticleList[pp]][0];
			s[1]+=ParticlePosition[ParticleList[pp]][1];
			s[2]+=ParticlePosition[ParticleList[pp]][2];
			NextnodeFromParticle[pp]=sib;
		}
		else
		{
			if(len>=HBTConfig.TreeNodeResolution) 
				UpdateInternalNodes(pp,sib,0.5*len);//only divide if above resolution; otherwise 
			//we didn't divide the node seriouly so we don't have finer node length
			else
				UpdateInternalNodes(pp,sib,len);
			mass+=Nodes[pp].way.mass;
			s[0]+=Nodes[pp].way.s[0]*Nodes[pp].way.mass;
			s[1]+=Nodes[pp].way.s[1]*Nodes[pp].way.mass;
			s[2]+=Nodes[pp].way.s[2]*Nodes[pp].way.mass;
		}
		Nodes[no].way.mass=mass;
		Nodes[no].way.s[0]=s[0]/mass;
		Nodes[no].way.s[1]=s[1]/mass;
		Nodes[no].way.s[2]=s[2]/mass;
}

HBTInt OctTree_t::Build(HBTInt num_part, HBTInt* particles, HBTReal part_pos[][3])
{/* build tree for a list of particles
  * particles[] contain the list of particle indices, 
  * part_pos[][3] is the position for all the particles
  * so that the position of particles[i] can be obtained as part_pos[particles[i]][]*/
	HBTInt NumNids,numnodes;
	HBTInt sub,subid,i,j,nodeid;
	double center[3], lenhalf;
	double xmin[3], xmax[3],Center[3], Len,Lenhalf;

	NumberOfParticles=num_part;
	ParticlePosition=part_pos;
	ParticleList=particles;

	Allocate(num_part);
	
	/* find enclosing rectangle */
  for(j = 0; j < 3; j++)
    xmin[j] = xmax[j] = ParticlePosition[ParticleList[0]][j];

  for(i = 1; i < NumberOfParticles; i++)
    for(j = 0; j < 3; j++)
      {
	if(ParticlePosition[ParticleList[i]][j] > xmax[j])
	  xmax[j] = ParticlePosition[ParticleList[i]][j];
	else if(ParticlePosition[ParticleList[i]][j] < xmin[j])
	  xmin[j] = ParticlePosition[ParticleList[i]][j];
      }

  /* determine maxmimum extension */
  for(j = 1, Len = xmax[0] - xmin[0]; j < 3; j++)
    if((xmax[j] - xmin[j]) > Len)
      Len = xmax[j] - xmin[j];

  for(j = 0; j < 3; j++)
    Center[j] = 0.5 * (xmax[j] + xmin[j]);
  
  Lenhalf=0.5*Len;
MaxNumberOfNodes=MaxNumberOfCells+NumberOfParticles;

Nodes= Cells-NumberOfParticles;	/* select first node */


nodeid = NumberOfParticles;	/* id used to distinguish whether it's internal node or particle*/
NumNids=NumberOfParticles+1;	
	/* create an empty  root node  */
  for(i = 0; i < 8; i++)
Cells->sons[i] = -1;

  for(i = 0; i < NumberOfParticles; i++)	/* insert all  particles */
	{
	  nodeid = NumberOfParticles ;	/* select index of first node in tree */
	    lenhalf = Lenhalf;
	  for(j = 0; j < 3; j++)
	    center[j] = Center[j];

	  while(1)
		{
			  //len = lenhalf;
			//fprintf(logfile,"%f\n",len);
			  lenhalf *= 0.5;//halflen for the to-be-found subnode
			  sub = 0;
	      if(ParticlePosition[ParticleList[i]][0] > center[0])
		{
		  center[0] += lenhalf;//subcenter
		  sub += 1;//sub index
		}
	      else
		{
		  center[0] -= lenhalf;
		}
	      if(ParticlePosition[ParticleList[i]][1] > center[1])
		{
		  center[1] += lenhalf;
		  sub += 2;
		}
	      else
		{
		  center[1] -= lenhalf;
		}
	      if(ParticlePosition[ParticleList[i]][2] > center[2])
		{
		  center[2] += lenhalf;
		  sub += 4;
		}
	      else
		{
		  center[2] -= lenhalf;
		}
		
		subid=Nodes[nodeid].sons[sub];
		if(subid<0)//an empty node, insert particle as leaf
			{
			Nodes[nodeid].sons[sub]=i;
			break;//finished for this particle, begin to insert a new particle
			}
		else if(subid<NumberOfParticles)//a particle node, upgrade the node to internal
			{
			Nodes[nodeid].sons[sub]=NumNids;//create a new node;
			nodeid=NumNids;//take over the new nodeid
			NumNids++;
			if(NumNids >= MaxNumberOfNodes)
			{
			  
			  fprintf(stderr,"maximum number %zd of tree-nodes reached.\n", MaxNumberOfCells);
			  fprintf(stderr,"for particle "HBTIFMT"\n", i);
			  exit(1);
			}
			for(sub=0;sub<8;sub++)//initialize new node
				Nodes[nodeid].sons[sub]=-1;
			/*insert that subid into this new node*/
			//what if the two particles are too near? 
			//unnecessary to divide too fine, just get rid of one by random insertion. 
			 if(lenhalf < HBTConfig.TreeNodeResolutionHalf)
				{
				/* seems like we're dealing with particles   
				* at identical locations. randomize 
				* sub index (well below gravitational softening scale).
				* to introduce some ambiguity so that we don't get jammed!*/
				sub = (HBTInt) (8.0 * drand48());
				if(sub >= 8)
				sub = 7;
				//~ fprintf(logfile,"len=%g Len=%g sub=%d  i=%d (%g|%g|%g)\n",
					//~ lenhalf*2, Len, sub, i, ParticlePosition[ParticleList[i]][0], ParticlePosition[ParticleList[i]][1], ParticlePosition[ParticleList[i]][2]);
				}
			 else
				{
				sub=0;
				if(ParticlePosition[ParticleList[subid]][0] > center[0])
					sub += 1;
				if(ParticlePosition[ParticleList[subid]][1] > center[1])
					sub += 2;
				if(ParticlePosition[ParticleList[subid]][2] > center[2])
					sub += 4;
				}	
			Nodes[nodeid].sons[sub]=subid;//the disturbing particle inserted
			}
		else nodeid=subid;//an internal node,take over it;
		}
	}
	
	numnodes=NumNids-NumberOfParticles;
	//~ FilledFraction=(HBTReal)numnodes/MaxNumberOfCells;
	//~ #pragma omp critical
	//~ if(FilledFraction>MaxNodeFilledFraction) MaxNodeFilledFraction=FilledFraction;
	//~ fprintf(logfile,"used %d nodes out of allocated %d. (filled fraction %g)\n",
	 //~ numnodes, MaxNumberOfCells, (double)numnodes / MaxNumberOfCells);	
	/* finished inserting, now update for walk*/
	UpdateInternalNodes(NumberOfParticles , -1, Len);/*insert sibling and next infomation*/
	
	return numnodes; 
}

double OctTree_t::EvaluatePotential(HBTReal targetPos[3])
{
  OctTreeCell_t *nop = 0;
  HBTInt no;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pot, pos_x, pos_y, pos_z;

  pos_x = targetPos[0];
  pos_y = targetPos[1];
  pos_z = targetPos[2];

  h = 2.8 * HBTConfig.SofteningHalo;
  h_inv = 1.0 / h;

  pot = 0;

  no = NumberOfParticles;//start from root node

  while(no >= 0)
    {
      if(no < NumberOfParticles)		/* single particle */
	{
	  dx = ParticlePosition[ParticleList[no]][0] - pos_x;
	  dy = ParticlePosition[ParticleList[no]][1] - pos_y;
	  dz = ParticlePosition[ParticleList[no]][2] - pos_z;
	  if(HBTConfig.PeriodicBoundaryOn)
	  {
	  dx=NEAREST(dx);
	  dy=NEAREST(dy);
	  dz=NEAREST(dz);
	  }
	  mass = 1;
	  no = NextnodeFromParticle[no];
	      r2 = dx * dx + dy * dy + dz * dz;	
	}
      else
	{
	  nop = &Nodes[no];
	  dx = nop->way.s[0] - pos_x;
	  dy = nop->way.s[1] - pos_y;
	  dz = nop->way.s[2] - pos_z;
	  if(HBTConfig.PeriodicBoundaryOn)
	  {
	  dx=NEAREST(dx);
	  dy=NEAREST(dy);
	  dz=NEAREST(dz);
	  }
	  mass = nop->way.mass;
	  r2 = dx * dx + dy * dy + dz * dz;
		/* we have an internal node. Need to check opening criterion */
	  if((nop->way.len * nop->way.len )>( r2 * HBTConfig.TreeNodeOpenAngleSquare))
	    {
	      /* open cell */
	      no = nop->way.nextnode;
	      continue;
	    }
	  no = nop->way.sibling;	/* node can be used */
	}
 
      r = sqrt(r2);

      if(r >= h)
	pot -= mass / r;
      else
	{
	  u = r * h_inv;

	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp =
	      -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						   u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

	  pot += mass * h_inv * wp;
	}
    }

  return pot;
}

#define MIN_NUM_CELLS_OF_TREE 500
#define TreeAllocFactor 1. /* a value of 2 should be more than sufficient*/
void OctTree_t::Allocate(size_t num_part)
{
  MaxNumberOfCells =TreeAllocFactor*num_part;
  if(MaxNumberOfCells<MIN_NUM_CELLS_OF_TREE) MaxNumberOfCells=MIN_NUM_CELLS_OF_TREE;
  Cells=new OctTreeCell_t[MaxNumberOfCells+1];

  NextnodeFromParticle=new HBTInt[num_part];
}
void OctTree_t::Clear()
{
  if(NumberOfParticles)
  {
	delete [] NextnodeFromParticle;
	NumberOfParticles=0;
  }
  if(MaxNumberOfCells)
  {
	delete [] Cells;
	MaxNumberOfCells=0;
  }
  MaxNumberOfNodes=0;
}













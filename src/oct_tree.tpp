// #include <cstdlib>
// #include <cstdio>
// #include <cmath>
// #include <cstring>

#include "mymath.h"
#include "config_parser.h"

template <class T>
HBTInt OctTree_t<T>::Build(const Snapshot_t &snapshot, HBTInt num_part)
/* build tree for a snapshot (or SnapshotView); automatically resize memory if necessary.
  * if num_part>0 is given, then only use the first num_part particles in the snapshot
  */
{
	HBTInt NumNids,numnodes;
	HBTInt sub,subid,i,j,nodeid;
	double center[3], lenhalf;
	double xmin[3], xmax[3],Center[3], Len,Lenhalf;

	if(!num_part) num_part=snapshot.size();
	if(num_part>MaxNumberOfParticles)
	{
	  Clear();
	  Reserve(num_part);
	}
	
	NumberOfParticles=num_part;
	Snapshot=&snapshot;
	
	/* find enclosing rectangle */
  for(j = 0; j < 3; j++)
    xmin[j] = xmax[j] = Snapshot->GetComovingPosition(0)[j];

  for(i = 1; i < NumberOfParticles; i++)
    for(j = 0; j < 3; j++)
      {
	if(Snapshot->GetComovingPosition(i)[j] > xmax[j])
	  xmax[j] = Snapshot->GetComovingPosition(i)[j];
	else if(Snapshot->GetComovingPosition(i)[j] < xmin[j])
	  xmin[j] = Snapshot->GetComovingPosition(i)[j];
      }

  /* determine maxmimum extension */
  for(j = 1, Len = xmax[0] - xmin[0]; j < 3; j++)
    if((xmax[j] - xmin[j]) > Len)
      Len = xmax[j] - xmin[j];

  for(j = 0; j < 3; j++)
    Center[j] = 0.5 * (xmax[j] + xmin[j]);
  
  Lenhalf=0.5*Len;
MaxNodeId=MaxNumberOfCells+NumberOfParticles;

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
	      if(Snapshot->GetComovingPosition(i)[0] > center[0])
		{
		  center[0] += lenhalf;//subcenter
		  sub += 1;//sub index
		}
	      else
		{
		  center[0] -= lenhalf;
		}
	      if(Snapshot->GetComovingPosition(i)[1] > center[1])
		{
		  center[1] += lenhalf;
		  sub += 2;
		}
	      else
		{
		  center[1] -= lenhalf;
		}
	      if(Snapshot->GetComovingPosition(i)[2] > center[2])
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
			if(NumNids >= MaxNodeId)
			{
			  cerr<<"maximum number "<<MaxNumberOfCells<<" of tree-nodes reached for particle "<<i<<endl;
			  exit(1);
			}
			for(sub=0;sub<8;sub++)//initialize new node
				Nodes[nodeid].sons[sub]=-1;
// 			copyXYZ(Nodes[nodeid].way.s, center);//DO NOT DO IT HERE: corrupting the union data
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
					//~ lenhalf*2, Len, sub, i, Snapshot->GetComovingPosition(i][0], Snapshot->GetComovingPosition(i][1], Snapshot->GetComovingPosition(i][2]);
				}
			 else
				{
				sub=0;
				if(Snapshot->GetComovingPosition(subid)[0] > center[0])
					sub += 1;
				if(Snapshot->GetComovingPosition(subid)[1] > center[1])
					sub += 2;
				if(Snapshot->GetComovingPosition(subid)[2] > center[2])
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
	UpdateInternalNodes(NumberOfParticles , -1, Len, Center);/*insert sibling and next infomation*/
	
	return numnodes; 
}

template <class T>
void OctTree_t<T>::Reserve(const size_t max_num_part)
/* allocate tree memory to hold a maximum of max_num_part particles */
{
  MaxNumberOfParticles=max_num_part;
  MaxNumberOfCells =HBTConfig.TreeAllocFactor*MaxNumberOfParticles;
  if(MaxNumberOfCells<HBTConfig.TreeMinNumOfCells) MaxNumberOfCells=HBTConfig.TreeMinNumOfCells;
  Cells=new OctTreeCell_t[MaxNumberOfCells+1];

  NextnodeFromParticle=new HBTInt[MaxNumberOfParticles];
}

template <class T>
void OctTree_t<T>::Clear()
{
  if(MaxNumberOfParticles)
  {
	delete [] NextnodeFromParticle;
	MaxNumberOfParticles=0;
	NumberOfParticles=0;
  }
  if(MaxNumberOfCells)
  {
	delete [] Cells;
	MaxNumberOfCells=0;
	MaxNodeId=0;
  }
}










#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "mymath.h"
#include "config_parser.h"
#include "gravity_tree.h"

template <class T>
inline void VectorAdd(double x[3], const T &y, double weight)
{
  x[0]+=y[0]*weight;
  x[1]+=y[1]*weight;
  x[2]+=y[2]*weight;
}

inline void GravityTree_t::FillNodeCenter(HBTInt nodeid, const double center[3], double CoM[3], double mass)
{
    Nodes[nodeid].way.s[0]=CoM[0]/mass;
    Nodes[nodeid].way.s[1]=CoM[1]/mass;
    Nodes[nodeid].way.s[2]=CoM[2]/mass;
}

void GravityTree_t::ProcessNode(HBTInt nodeid, HBTInt nextid, int sonid, double &mass, double CoM[3], double len, const double center[3])
{
  if(nodeid<NumberOfParticles)
      {
	double thismass=Snapshot->GetMass(nodeid);
	mass+=thismass;
	VectorAdd(CoM, Snapshot->GetComovingPosition(nodeid), thismass);
	
	NextnodeFromParticle[nodeid]=nextid;
      }
      else
      {
	double thismass=Nodes[nodeid].way.mass;
	mass+=thismass;
	VectorAdd(CoM, Nodes[nodeid].way.s, thismass);
	
	if(len>=HBTConfig.TreeNodeResolution)//only divide if above resolution;
	  UpdateInternalNodes(nodeid, nextid, len/2., center); 
	else
	  UpdateInternalNodes(nodeid, nextid, len, center);//otherwise we don't divide the node seriouly so we don't have finer node length
      }
}

void GravityTree_t::UpdateInternalNodes(HBTInt no, HBTInt sib, double len, const double center[3])
{
  HBTInt p,pp,sons[8];
  int j,jj,i;
  double mass=0., thismass;
  double CoM[3]={0.};
  
  for(j=0;j<8;j++)
    sons[j]=Nodes[no].sons[j];//backup sons
  Nodes[no].way.len=len;
  Nodes[no].way.sibling=sib;
  for(i=0;sons[i]<0;i++);//find first son
  jj=i;
  pp=sons[jj];
  Nodes[no].way.nextnode=pp;
  for(i++;i<8;i++)//find sons in pairs,ie. find sibling
  {
    if(sons[i]>=0)//ok, found a sibling
    {
      j=jj;
      p=pp;
      jj=i;
      pp=sons[jj];
      ProcessNode(p, pp, j, mass, CoM, len, center);
    }
  }
  ProcessNode(pp, sib, jj, mass, CoM, len, center);
  Nodes[no].way.mass=mass;
  FillNodeCenter(no, center, CoM, mass);
}

double GravityTree_t::EvaluatePotential(const HBTxyz &targetPos, const HBTReal targetMass)
/*return specific physical potential, GM/Rphysical. 
 * targetPos[] is comoving.
 * if targetMass!=0, then the self-potential from targetMass is excluded. 
 * do not set targetMass (i.e., keep to 0.) if target is outside the particlelist of tree*/
{
  bool IsPeriodic=HBTConfig.PeriodicBoundaryOn;
  OctTreeCell_t *nop = 0;
  HBTInt no;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pot, pos_x, pos_y, pos_z;

  pos_x = targetPos[0];
  pos_y = targetPos[1];
  pos_z = targetPos[2];

  h = 2.8 * HBTConfig.SofteningHalo;
  h_inv = 1.0 / h;

  pot=targetMass/HBTConfig.SofteningHalo; //to cancle out the self-potential added during tree walk.

  no = NumberOfParticles;//start from root node

  while(no >= 0)
    {
      if(no < NumberOfParticles)		/* single particle */
	{
	  auto &pos=Snapshot->GetComovingPosition(no);
	  dx = pos[0] - pos_x;
	  dy = pos[1] - pos_y;
	  dz = pos[2] - pos_z;
	  if(IsPeriodic)
	  {
	  dx=NEAREST(dx);
	  dy=NEAREST(dy);
	  dz=NEAREST(dz);
	  }
	  mass = Snapshot->GetMass(no);
	  no = NextnodeFromParticle[no];
	      r2 = dx * dx + dy * dy + dz * dz;	
	}
      else
	{
	  nop = &Nodes[no];
	  dx = nop->way.s[0] - pos_x;
	  dy = nop->way.s[1] - pos_y;
	  dz = nop->way.s[2] - pos_z;
	  if(IsPeriodic)
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
    
  return pot*PhysicalConst::G/Snapshot->Cosmology.ScaleFactor;
}

double GravityTree_t::BindingEnergy(const HBTxyz& targetPos, const HBTxyz& targetVel, const HBTxyz& refPos, const HBTxyz& refVel, const HBTReal targetMass)
/* return specific binding energy 
 * input Pos comoving, Vel physical
 * targetMass optional, can be set to exclude self-potential if target is contained in the tree*/
{
	  double pot=EvaluatePotential(targetPos, targetMass);
	  HBTxyz dv;
	  Snapshot->RelativeVelocity(targetPos, targetVel, refPos, refVel, dv);
	  return VecNorm(dv)*0.5+pot;
}

template class OctTree_t<GravityTreeCell_t>;//to wake up the functions for this type; trick!

#ifdef TEST_gravity_tree
#include "snapshot.h"
#include "halo.h"

int main(int argc, char **argv)
{
  HBTConfig.ParseConfigFile("../configs/AqA5.conf");
  HBTInt isnap=HBTConfig.MinSnapshotIndex;
  ParticleSnapshot_t snapshot;
  snapshot.Load(isnap);
  
  HaloSnapshot_t halo;
  halo.Load(isnap);
  halo.ParticleIdToIndex(snapshot);

  halo.AverageCoordinates();
  Halo_t::ParticleList_t &P=halo.Halos[0].Particles;
  
  SnapshotView_t treesnap(P, snapshot);
  
  OctTree_t tree;
  tree.Reserve(2);
  tree.Build(treesnap);
  
  for(HBTInt i=0;i<P.size();i++)
  {
	double E=tree.BindingEnergy(snapshot.GetComovingPosition(P[i]), snapshot.GetPhysicalVelocity(P[i]),
								halo.Halos[0].ComovingPosition, halo.Halos[0].PhysicalVelocity, snapshot.GetMass(P[i]));
	cout<<i<<" : "<<E<<endl;
  }
  
  return 0;
}
#endif










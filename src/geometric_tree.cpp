#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "mymath.h"
#include "config_parser.h"
#include "geometric_tree.h"

inline void shift_center(const double oldcenter[3], int son, double delta, double newcenter[3])
{
  for(int dim=0;dim<3;dim++)
  {
    int bit=get_bit(son, dim);
    if(bit)
      newcenter[dim]=oldcenter[dim]+delta;
    else
      newcenter[dim]=oldcenter[dim]-delta;
  }
}

void GeoTree_t::ProcessNode(HBTInt nodeid, HBTInt nextid, int sonid, HBTInt &mass, double len, const double center[3])
{     
      if(nodeid<NumberOfParticles)
      {
	mass++;
	NextnodeFromParticle[nodeid]=nextid;
      }
      else
      {	
	if(len>=HBTConfig.TreeNodeResolution)//only divide if above resolution;
	{
	  double newcenter[3];
	  shift_center(center, sonid, len/4., newcenter);
	  UpdateInternalNodes(nodeid, nextid, len/2., newcenter);
	}
	else
	  UpdateInternalNodes(nodeid, nextid, len, center);//otherwise we don't divide the node seriouly so we don't have finer node length

	mass+=Nodes[nodeid].way.mass;//update mass after updating internal nodes
      }
}

inline void GeoTree_t::FillNodeCenter(HBTInt nodeid, const double center[3])
{
  copyXYZ(Nodes[nodeid].way.s, center);
}

void GeoTree_t::UpdateInternalNodes(HBTInt no, HBTInt sib, double len, const double center[3])
{
  HBTInt p,pp,sons[8];
  int j,jj,i;
  HBTInt mass=0;
  
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
      ProcessNode(p, pp, j, mass, len, center);
    }
  }
  ProcessNode(pp, sib, jj, mass, len, center);
  Nodes[no].way.mass=mass;
  FillNodeCenter(no, center);
}

inline HBTReal GuessNeighbourRange(HBTInt n_neighbours, HBTReal number_density_guess)
{
  return pow(3 * n_neighbours / (4 * 3.141593) / number_density_guess, 1.0 / 3);
}

HBTInt GeoTree_t::NearestNeighbour(const HBTxyz & cen, HBTReal rguess)
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

void GeoTree_t::Search(const HBTxyz & searchcenter, HBTReal radius, vector <LocatedParticle_t> &founds)
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
	  if(dx > radius || dx < -radius)
	    continue;
	  
	  double dy = pos[1] - y0;
	  if(IsPeriodic) dy=NEAREST(dy);
	  if(dy > radius || dy < -radius)
	    continue;

	  double dz = pos[2] - z0;
	  if(IsPeriodic) dz=NEAREST(dz);
	  if(dz > radius || dz < -radius)
	    continue;

	  double r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	      founds.emplace_back(Snapshot->GetMemberId(pid), r2);
	}
      else
	{
	  auto &node = Nodes[node_id];

	  node_id = node.way.sibling;	/* in case the node can be discarded */
	  double rmax=node.way.len/2.;
	  rmax+=radius;
	  
	  auto &pos=node.way.s;
	  double dx = pos[0] - x0;
	  if(IsPeriodic) dx=NEAREST(dx);
	  if(dx > rmax || dx < -rmax)
	    continue;
	  
	  double dy = pos[1] - y0;
	  if(IsPeriodic) dy=NEAREST(dy);
	  if(dy > rmax || dy < -rmax)
	    continue;

	  double dz = pos[2] - z0;
	  if(IsPeriodic) dz=NEAREST(dz);
	  if(dz > rmax || dz < -rmax)
	    continue;

	  node_id = node.way.nextnode;	/* ok, we need to open the node */
	}
    }
}

#define SPH_DENS_NGB 64
double GeoTree_t::SphDensity(const HBTxyz &cen, HBTReal & hguess)
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

template class OctTree_t<GeoTreeCell_t>;//to wake up the functions for this type; trick!

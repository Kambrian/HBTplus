#include <iostream>
#include <new>
#include <algorithm>
#include <unordered_set>
#include <omp.h>

#include "datatypes.h"
#include "snapshot_number.h"
#include "subhalo.h"

#define NumPartCoreMax 20
#define DeltaCrit 2.

struct SubHelper_t
{
  HBTInt HostTrackId;
  bool IsMerged;
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  float ComovingSigmaR;
  float PhysicalSigmaV;
  void BuildPosition(const Subhalo_t &sub);
  void BuildVelocity(const Subhalo_t &sub);
  SubHelper_t(): HostTrackId(-1), IsMerged(false)
  {
  }
};

void SubHelper_t::BuildPosition(const Subhalo_t &sub)
{
  if(0==sub.Nbound)
  {
    ComovingSigmaR=0.;
    return;
  }
  if(1==sub.Nbound) 
  {
    ComovingSigmaR=0.;
    copyHBTxyz(ComovingPosition, sub.Particles[0].ComovingPosition);
    return;
  }
  
  HBTInt NumPart=sub.Nbound;
  if(NumPart>NumPartCoreMax) NumPart=NumPartCoreMax;
  
  HBTInt i,j;
  double sx[3],sx2[3], origin[3],msum;
    
  sx[0]=sx[1]=sx[2]=0.;
  sx2[0]=sx2[1]=sx2[2]=0.;
  msum=0.;
  if(HBTConfig.PeriodicBoundaryOn)
    for(j=0;j<3;j++)
      origin[j]=sub.Particles[0].ComovingPosition[j];
    
    for(i=0;i<NumPart;i++)
    {
      HBTReal m=sub.Particles[i].Mass;
      msum+=m;
      for(j=0;j<3;j++)
      {
	double dx;
	if(HBTConfig.PeriodicBoundaryOn)
	  dx=NEAREST(sub.Particles[i].ComovingPosition[j]-origin[j]);
	else
	  dx=sub.Particles[i].ComovingPosition[j];
	sx[j]+=dx*m;
	sx2[j]+=dx*dx*m;
      }
    }
    
    for(j=0;j<3;j++)
    {
      sx[j]/=msum;
      sx2[j]/=msum;
      ComovingPosition[j]=sx[j];
      if(HBTConfig.PeriodicBoundaryOn) ComovingPosition[j]+=origin[j];
      sx2[j]-=sx[j]*sx[j];
    }
    ComovingSigmaR=sqrt(sx2[0]+sx2[1]+sx2[2]);
}
void SubHelper_t::BuildVelocity(const Subhalo_t &sub)
{
  if(0==sub.Nbound)
  {
    PhysicalSigmaV=0.;
    return;
  }
  if(1==sub.Nbound) 
  {
    PhysicalSigmaV=0.;
    copyHBTxyz(PhysicalVelocity, sub.Particles[0].PhysicalVelocity);
    return;
  }
  
  HBTInt NumPart=sub.Nbound;
  if(NumPart>NumPartCoreMax) NumPart=NumPartCoreMax;
  
  HBTInt i,j;
  double sx[3],sx2[3],msum;
  
  sx[0]=sx[1]=sx[2]=0.;
  sx2[0]=sx2[1]=sx2[2]=0.;
  msum=0.;
  
  for(i=0;i<NumPart;i++)
  {
    HBTReal m=sub.Particles[i].Mass;
    msum+=m;
    for(j=0;j<3;j++)
    {
      double dx;
      dx=sub.Particles[i].PhysicalVelocity[j];
      sx[j]+=dx*m;
      sx2[j]+=dx*dx*m;
    }
  }
  
  for(j=0;j<3;j++)
  {
    sx[j]/=msum;
    sx2[j]/=msum;
    PhysicalVelocity[j]=sx[j];
    sx2[j]-=sx[j]*sx[j];
  }
  PhysicalSigmaV=sqrt(sx2[0]+sx2[1]+sx2[2]);
}

float SinkDistance(const Subhalo_t &sat, const SubHelper_t &cen)
{
  float d=PeriodicDistance(cen.ComovingPosition, sat.ComovingMostBoundPosition);
  float v=Distance(cen.PhysicalVelocity, sat.PhysicalMostBoundVelocity);
  return d/cen.ComovingSigmaR+v/cen.PhysicalSigmaV;
}

void DetectTraps(vector <Subhalo_t> &Subhalos, vector <SubHelper_t> &Helpers)
{  
	#pragma omp  for schedule(dynamic,1)
	for(HBTInt i=0;i<Subhalos.size();i++)
	{
// 	  if(Subhalos[i].Nbound<=1) continue;//skip orphans? no.
	  if(Subhalos[i].IsTrapped()) continue;
	  HBTInt HostId=Helpers[i].HostTrackId;
	  while(HostId>=0)
	  {
	    if(Subhalos[HostId].Nbound>1)//avoid orphans or nulls as hosts
	    {
	      float delta=SinkDistance(Subhalos[i], Helpers[HostId]);
	      if(delta<DeltaCrit)
	      {
		Subhalos[i].SinkTrackId=HostId;//these are local ids for the merging tracks. Those already merged ones retain their global ids.
		if(Subhalos[i].Nbound>1) //only need to unbind if a real sub sinks
		  Helpers[HostId].IsMerged=true;
		break;
	      }
	    }
	    HostId=Helpers[HostId].HostTrackId;
	  }
	}
}

void FillHostTrackIds(vector <SubHelper_t> &Helpers, const vector <Subhalo_t> &Subhalos)
{
  #pragma omp for schedule(dynamic,1)
  for(HBTInt i=0;i<Subhalos.size();i++) 
  {
    auto &nest=Subhalos[i].NestedSubhalos;
    for(auto &&subid: nest)
      Helpers[subid].HostTrackId=i;
  }
}
void FillCores(vector <SubHelper_t> &Helpers, const vector <Subhalo_t> &Subhalos)
{
  #pragma omp for schedule(dynamic,1) 
    for(HBTInt i=0;i<Subhalos.size();i++)
    {
      Helpers[i].BuildPosition(Subhalos[i]);
      Helpers[i].BuildVelocity(Subhalos[i]);
    }
}
void FillHelpers(vector <SubHelper_t> &Helpers, const vector <Subhalo_t> &Subhalos)
{
  FillHostTrackIds(Helpers, Subhalos);
  FillCores(Helpers, Subhalos);
}
void SubhaloSnapshot_t::MergeSubhalos()
{
  HBTInt NumHalos=MemberTable.SubGroups.size();
  vector <SubHelper_t> Helpers(Subhalos.size());
  
  #pragma omp parallel
  {
    GlueHeadNests();
    FillHelpers(Helpers, Subhalos);
    
    DetectTraps(Subhalos, Helpers);
  }
  
  if(HBTConfig.MergeTrappedSubhalos)
  {
    #pragma omp parallel for schedule(dynamic,1)
    for(HBTInt grpid=0;grpid<NumHalos;grpid++)
      if(MemberTable.SubGroups[grpid].size())
	MergeRecursive(MemberTable.SubGroups[grpid][0]);
    #pragma omp parallel for schedule(dynamic,1) if(ParallelizeHaloes)
    for(HBTInt subid=0;subid<Subhalos.size();subid++)
      if(Helpers[subid].IsMerged)
	Subhalos[subid].Unbind(*this);
    #pragma omp parallel for
    for(HBTInt subid=0;subid<Subhalos.size();subid++)
      if(Helpers[subid].IsMerged)
	Subhalos[subid].TruncateSource();
  }
  
  Helpers.clear();
  #pragma omp parallel
  UnglueHeadNests();
}

void SubhaloSnapshot_t::MergeRecursive(HBTInt subid)
{
  auto &sat=Subhalos[subid];
  for(auto &&nestid: sat.NestedSubhalos)
    MergeRecursive(nestid);
  if(sat.IsTrapped()&&sat.IsAlive())
  {
    sat.MergeTo(Subhalos[sat.SinkTrackId]);
    sat.SnapshotIndexOfDeath=GetSnapshotIndex();
  }
}


/*
struct ParticleHasher_t //to be passed as a template parameter to unordered_set
{
  size_t operator()(const Particle_t& p) const
  {
    return p.Id;
  }
};
*/
void Subhalo_t::MergeTo(Subhalo_t &host)
{
  if(Nbound<=1) return; //skip orphans and nulls
  
  #ifndef INCLUSIVE_MASS
  HBTInt np_max=host.Particles.size()+Particles.size();
  unordered_set <HBTInt> UniqueIds(np_max);
  for(auto &&p: host.Particles)   UniqueIds.insert(p.Id);
  host.Particles.reserve(np_max);
  for(auto &&p: Particles)
    if(UniqueIds.insert(p.Id).second)//inserted, meaning not excluded
      host.Particles.push_back(p);
  host.Nbound+=Nbound;
  #endif
  
  Particles.resize(1);
  Nbound=1;
}

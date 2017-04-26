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

struct SubCore_t
{
  HBTxyz ComovingPosition;
  HBTxyz PhysicalVelocity;
  float ComovingSigmaR;
  float PhysicalSigmaV;
  void BuildPosition(const Subhalo_t &sub, const Snapshot_t &snap);
  void BuildVelocity(const Subhalo_t &sub, const Snapshot_t &snap);
};

void SubCore_t::BuildPosition(const Subhalo_t &sub, const Snapshot_t &snap)
{
  if(0==sub.Nbound)
  {
    ComovingSigmaR=0.;
    return;
  }
  if(1==sub.Nbound) 
  {
    ComovingSigmaR=0.;
    copyHBTxyz(ComovingPosition, snap.GetComovingPosition(sub.Particles[0]));
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
      origin[j]=snap.GetComovingPosition(sub.Particles[0])[j];
    
    for(i=0;i<NumPart;i++)
    {
      HBTReal m=snap.GetMass(sub.Particles[i]);
      msum+=m;
      for(j=0;j<3;j++)
      {
	double dx;
	if(HBTConfig.PeriodicBoundaryOn)
	  dx=NEAREST(snap.GetComovingPosition(sub.Particles[i])[j]-origin[j]);
	else
	  dx=snap.GetComovingPosition(sub.Particles[i])[j];
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
void SubCore_t::BuildVelocity(const Subhalo_t &sub, const Snapshot_t &snap)
{
  if(0==sub.Nbound)
  {
    PhysicalSigmaV=0.;
    return;
  }
  if(1==sub.Nbound) 
  {
    PhysicalSigmaV=0.;
    copyHBTxyz(PhysicalVelocity, snap.GetPhysicalVelocity(sub.Particles[0]));
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
    HBTReal m=snap.GetMass(sub.Particles[i]);
    msum+=m;
    for(j=0;j<3;j++)
    {
      double dx;
      dx=snap.GetPhysicalVelocity(sub.Particles[i])[j];
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

float SinkDistance(const Subhalo_t &sat, const SubCore_t &cen)
{
  float d=PeriodicDistance(cen.ComovingPosition, sat.ComovingMostBoundPosition);
  float v=Distance(cen.PhysicalVelocity, sat.PhysicalMostBoundVelocity);
  return d/cen.ComovingSigmaR+v/cen.PhysicalSigmaV;
}

void DetectTraps(vector <Subhalo_t> &Subhalos, const vector <SubCore_t> &Cores, int isnap)
{  
	#pragma omp  for
	for(HBTInt i=0;i<Subhalos.size();i++)
	{
// 	  if(Subhalos[i].Nbound<=1) continue;//skip orphans
	  if(Subhalos[i].IsTrapped)
	  {
// 	    Subhalos[i].Delta=SinkDistance(Subhalos[i], Subhalos[Subhalos[i].SinkTrackId]);
	    continue;
	  }
	  HBTInt HostId=Subhalos[i].HostTrackId;
	  float deltaMin=-1.;
	  HBTInt SinkId=HostId;
	  while(HostId>=0)
	  {
	    if(Subhalos[HostId].Nbound>1)//skip orphans and nulls
	    {
	      float delta=SinkDistance(Subhalos[i], Cores[HostId]);
	      if(delta<deltaMin||deltaMin<0) 
	      {
		deltaMin=delta;
		SinkId=HostId;
	      } 
	      if(delta<DeltaCrit)
	      {
  // 	      Subhalos[i].Delta=delta;
		Subhalos[i].SinkTrackId=HostId;
  // 	      Subhalos[i].TrackDepthAtSink=Subhalos[i].Depth;
  // 	      Subhalos[i].SinkTrackDepthAtSink=Subhalos[HostId].Depth;
  // 	      Subhalos[i].DeltaAtSink=delta;
  // 	      Subhalos[i].MboundSink=Subhalos[i].Mbound;
  // 	      Subhalos[i].MratSink=Subhalos[i].Mbound/Subhalos[HostId].Mbound;
		Subhalos[i].IsTrapped=1;
		break;
	      }
	    }
	    HostId=Subhalos[HostId].HostTrackId;
	  }
	  if(!Subhalos[i].IsTrapped&&SinkId>=0) //bind to the minimum delta before sink
	  {
// 	    Subhalos[i].Delta=deltaMin;
	    Subhalos[i].SinkTrackId=SinkId;
	  }
	}
}

void SubhaloSnapshot_t::FillHostTrackIds()
{
      #pragma omp for
      for(HBTInt i=0;i<Subhalos.size();i++) 
      {
	Subhalos[i].HostTrackId=-1;
// 	Subhalos[i].Delta=-1.;//reinit
      }
      #pragma omp for
      for(HBTInt i=0;i<Subhalos.size();i++) 
      {
	auto &nest=Subhalos[i].NestedSubhalos;
	for(auto &&subid: nest)
	  Subhalos[subid].HostTrackId=i;
      }
}
void FillCores(vector <SubCore_t> &Cores, const vector <Subhalo_t> &Subhalos, const Snapshot_t &snap)
{
  #pragma omp for
    for(HBTInt i=0;i<Subhalos.size();i++)
    {
      Cores[i].BuildPosition(Subhalos[i], snap);
      Cores[i].BuildVelocity(Subhalos[i], snap);
    }
}
void SubhaloSnapshot_t::MergeSubhalos()
{
  HBTInt NumHalos=MemberTable.SubGroups.size();
  vector <int> membercount(NumHalos, 0);  
  vector <SubCore_t> Cores(Subhalos.size());
  
  #pragma omp parallel
  {
    #pragma omp for
    for(HBTInt haloid=0;haloid<NumHalos;haloid++)
    {//restore nest to the state during unbinding
	  auto &subgroup=MemberTable.SubGroups[haloid];
	  if(subgroup.size()==0) continue;
	  auto &nests=Subhalos[subgroup[0]].NestedSubhalos;
	  membercount[haloid]=nests.size();
	  auto &heads=MemberTable.SubGroupsOfHeads[haloid];
	  nests.insert(nests.end(), heads.begin()+1, heads.end());
    }
    FillHostTrackIds();
    
    FillCores(Cores, Subhalos, *SnapshotPointer);
    DetectTraps(Subhalos, Cores, GetSnapshotIndex());
    #pragma omp single
    Cores.clear();
    
    if(HBTConfig.MergeTrappedSubhalos)
    {
    #pragma omp for
    for(HBTInt subid=0;subid<Subhalos.size();subid++)
      Subhalos[subid].NumberOfMergers=0;
    #pragma omp for
    for(HBTInt grpid=0;grpid<NumHalos;grpid++)
      if(MemberTable.SubGroups[grpid].size())
	MergeRecursive(MemberTable.SubGroups[grpid][0]);
    #pragma omp for
    for(HBTInt subid=0;subid<Subhalos.size();subid++)
      if(Subhalos[subid].NumberOfMergers)
      {
	Subhalos[subid].Unbind(*SnapshotPointer);
	Subhalos[subid].TruncateSource();
      }
    }
    
    #pragma omp for
    for(HBTInt haloid=0;haloid<NumHalos;haloid++)
    {
	  auto &subgroup=MemberTable.SubGroups[haloid];
	  if(subgroup.size()) 
	    Subhalos[subgroup[0]].NestedSubhalos.resize(membercount[haloid]);//restore old satellite list
    }
  }
}

void SubhaloSnapshot_t::MergeRecursive(HBTInt subid)
{
  auto &sat=Subhalos[subid];
  for(auto &&nestid: sat.NestedSubhalos)
    MergeRecursive(nestid);
  if(sat.IsTrapped&&sat.SnapshotIndexOfDeath==SpecialConst::NullSnapshotId)
  {
    sat.MergeTo(Subhalos[sat.SinkTrackId]);
    if(sat.SnapshotIndexOfDeath==SpecialConst::NullSnapshotId)
      sat.SnapshotIndexOfDeath=GetSnapshotIndex();
  }
}

void Subhalo_t::MergeTo(Subhalo_t &host)
{
  if(Nbound<=1) return;//skip orphans and nulls
  host.NumberOfMergers+=NumberOfMergers+1;//including hierarchical (indirect) mergers.

  #ifndef INCLUSIVE_MASS
  unordered_set <HBTInt> UniqueParticles(host.Particles.begin(), host.Particles.end());
  UniqueParticles.insert(Particles.begin(), Particles.end());
  host.Particles.assign(UniqueParticles.begin(), UniqueParticles.end());
  host.Nbound+=Nbound;
  #endif
  
  Particles.resize(1);
  Nbound=1;
}

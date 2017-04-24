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

void Subhalo_t::CalcPositionCore(const Snapshot_t &snap)
{
  if(0==Nbound)
  {
    ComovingCoreSigmaR=0.;
    return;
  }
  if(1==Nbound) 
  {
    ComovingCoreSigmaR=0.;
    copyHBTxyz(ComovingCorePosition, snap.GetComovingPosition(Particles[0]));
    return;
  }
  
  HBTInt NumPart=Nbound;
  if(NumPart>NumPartCoreMax) NumPart=NumPartCoreMax;
  
  HBTInt i,j;
  double sx[3],sx2[3], origin[3],msum;
    
  sx[0]=sx[1]=sx[2]=0.;
  sx2[0]=sx2[1]=sx2[2]=0.;
  msum=0.;
  if(HBTConfig.PeriodicBoundaryOn)
    for(j=0;j<3;j++)
      origin[j]=snap.GetComovingPosition(Particles[0])[j];
    
    for(i=0;i<NumPart;i++)
    {
      HBTReal m=snap.GetParticleMass(Particles[i]);
      msum+=m;
      for(j=0;j<3;j++)
      {
	double dx;
	if(HBTConfig.PeriodicBoundaryOn)
	  dx=NEAREST(snap.GetComovingPosition(Particles[i])[j]-origin[j]);
	else
	  dx=snap.GetComovingPosition(Particles[i])[j];
	sx[j]+=dx*m;
	sx2[j]+=dx*dx*m;
      }
    }
    
    for(j=0;j<3;j++)
    {
      sx[j]/=msum;
      sx2[j]/=msum;
      ComovingCorePosition[j]=sx[j];
      if(HBTConfig.PeriodicBoundaryOn) ComovingCorePosition[j]+=origin[j];
      sx2[j]-=sx[j]*sx[j];
    }
    ComovingCoreSigmaR=sqrt(sx2[0]+sx2[1]+sx2[2]);
}
void Subhalo_t::CalcVelocityCore(const Snapshot_t &snap)
{
  if(0==Nbound)
  {
    PhysicalCoreSigmaV=0.;
    return;
  }
  if(1==Nbound) 
  {
    PhysicalCoreSigmaV=0.;
    copyHBTxyz(PhysicalCoreVelocity, snap.GetPhysicalVelocity(Particles[0]));
    return;
  }
  
  HBTInt NumPart=Nbound;
  if(NumPart>NumPartCoreMax) NumPart=NumPartCoreMax;
  
  HBTInt i,j;
  double sx[3],sx2[3],msum;
  
  sx[0]=sx[1]=sx[2]=0.;
  sx2[0]=sx2[1]=sx2[2]=0.;
  msum=0.;
  
  for(i=0;i<NumPart;i++)
  {
    HBTReal m=snap.GetParticleMass(Particles[i]);
    msum+=m;
    for(j=0;j<3;j++)
    {
      double dx;
      dx=snap.GetPhysicalVelocity(Particles[i])[j];
      sx[j]+=dx*m;
      sx2[j]+=dx*dx*m;
    }
  }
  
  for(j=0;j<3;j++)
  {
    sx[j]/=msum;
    sx2[j]/=msum;
    PhysicalCoreVelocity[j]=sx[j];
    sx2[j]-=sx[j]*sx[j];
  }
  PhysicalCoreSigmaV=sqrt(sx2[0]+sx2[1]+sx2[2]);
}

float SinkDistance(Subhalo_t &sat, Subhalo_t &cen)
{
  float d=PeriodicDistance(cen.ComovingCorePosition, sat.ComovingMostBoundPosition);
  float v=Distance(cen.PhysicalCoreVelocity, sat.PhysicalMostBoundVelocity);
  return d/cen.ComovingCoreSigmaR+v/cen.PhysicalCoreSigmaV;
}

void SubhaloSnapshot_t::DetectTraps()
{
  int isnap=GetSnapshotIndex();
  auto & Satellites=Subhalos;

  auto &SubGroups=MemberTable.SubGroups;
    
    #pragma omp parallel
    {
      //reinit host
      #pragma omp for
      for(HBTInt i=0;i<Satellites.size();i++) 
      {
	Satellites[i].HostTrackId=-1;
// 	Satellites[i].Delta=-1.;//reinit
	Satellites[i].CalcPositionCore(*SnapshotPointer);
	Satellites[i].CalcVelocityCore(*SnapshotPointer);
      }
      #pragma omp for
      for(HBTInt i=0;i<Satellites.size();i++) 
      {
	auto &nest=subsnap.Subhalos[i].NestedSubhalos;
	for(auto &&subid: nest)
	  Satellites[subid].HostTrackId=i;
      }
      	
	#pragma omp  for
	for(HBTInt i=0;i<Satellites.size();i++)
	{
// 	  if(Satellites[i].Nbound<=1) continue;//skip orphans
	  if(Satellites[i].HasSinked())
	  {
// 	    Satellites[i].Delta=SinkDistance(Subhalos[i], Satellites[Satellites[i].SinkTrackId]);
	    continue;
	  }
	  HBTInt HostId=Satellites[i].HostTrackId;
	  float deltaMin=-1.;
	  HBTInt SinkId=HostId;
	  while(HostId>=0)
	  {
	    float delta=SinkDistance(Subhalos[i], Satellites[HostId]);
	    if(delta<deltaMin||deltaMin<0) 
	    {
	      deltaMin=delta;
	      SinkId=HostId;
	    } 
	    if(delta<DeltaCrit)
	    {
// 	      Satellites[i].Delta=delta;
	      Satellites[i].SinkTrackId=HostId;
// 	      Satellites[i].TrackDepthAtSink=Satellites[i].TrackDepth;
// 	      Satellites[i].SinkTrackDepthAtSink=Satellites[HostId].TrackDepth;
// 	      Satellites[i].DeltaAtSink=delta;
// 	      Satellites[i].MboundSink=subsnap.Subhalos[i].Mbound;
// 	      Satellites[i].MratSink=subsnap.Subhalos[i].Mbound/subsnap.Subhalos[HostId].Mbound;
	      Satellites[i].SnapshotIndexOfSink=isnap;
	      break;
	    }
	    HostId=Satellites[HostId].HostTrackId;
	  }
	  if(!Satellites[i].HasSinked()&&SinkId>=0) //bind to the minimum delta before sink
	  {
// 	    Satellites[i].Delta=deltaMin;
	    Satellites[i].SinkTrackId=SinkId;
	  }
	}
    }
    
}

void SubhaloSnapshot_t::MergeSubhalos()
{
  int isnap=GetSnapshotIndex();
  for(HBTInt subid=0;subid<Subhalos.size();subid++)
  {
    auto &sat=Subhalos[subid];
    if((sat.Nbound>1)&&(sat.SnapshotIndexOfSink==isnap))
      Subhalos[sat.SinkTrackId].Merge(sat);
  }
}

Subhalo_t::Merge(Subhalo_t &sat)
{
  Particles.insert(Particles.begin()+Nbound, sat.Particles.begin(), sat.Particles.end());
  Nbound+=sat.Nbound;
  sat.Particles.resize(1);
  sat.Nbound=1;
}

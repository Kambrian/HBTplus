import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

MaxSnap=127
rootdir='/gpfs/data/jvbq85/HBT/data/AqA5/subcat2a/'

def NFWFunc(x):
  return np.log(1+x)-x/(1+x)

def NFWPot(sat, host, hubble, scale, XProxy):
  '''return potential and central potential'''
  Rs=host['RmaxComoving']/2.16
  c=host['R200CritComoving']/Rs
  R=distance(sat[XProxy], host[XProxy], axis=1)#input has to be arrays of halos
  x=R/Rs
  Pots=-100.*(hubble*Rs*scale)**2*c**3/NFWFunc(c)
  return Pots*np.log(1+x)/x, Pots

def getSub(trackId, isnap=MaxSnap):
  return h5py.File(rootdir+'/SubSnap_%03d.hdf5'%isnap,'r')['Subhalos'][trackId]

def getTrack(trackId):
  track=[];
  snaps=[]
  snapbirth=getSub(trackId)['SnapshotIndexOfBirth']
  for isnap in range(snapbirth, MaxSnap+1):
	  s=getSub(trackId, isnap)
	  track.append(s)
	  snaps.append(isnap)
  return np.array(track)

def getTrackPair(trackId):
  sattrack=[];
  hosttrack=[] #host at each snapshot
  Hosttrack=[] #final host track
  snaps=[]
  
  f=h5py.File(rootdir+'/SubSnap_%03d.hdf5'%MaxSnap,'r')
  sat=f['Subhalos'][trackId]
  HostTrackId=f['Membership/GroupedTrackIds'][sat['HostHaloId']][0]
  Host=f['Subhalos'][HostTrackId]
  snapbirth=max(sat['SnapshotIndexOfBirth'], Host['SnapshotIndexOfBirth'])
  #print snapbirth
  f.close()
  
  for isnap in range(snapbirth, MaxSnap+1):
	  f=h5py.File(rootdir+'/SubSnap_%03d.hdf5'%isnap,'r')
	  sat=f['Subhalos'][trackId]
	  hostTrackId=f['Membership/GroupedTrackIds'][sat['HostHaloId']][0]
	  host=f['Subhalos'][hostTrackId]
	  Host=f['Subhalos'][HostTrackId]
	  sattrack.append(sat)
	  hosttrack.append(host)
	  Hosttrack.append(Host)
	  snaps.append([isnap, f['HubbleParam'][...], f['ScaleFactor'][...]])
	  f.close()
	  
  return np.array(sattrack), np.array(Hosttrack), np.array(hosttrack), np.array(snaps)

def getHostTrackId(trackId):
  hosttrackId=[]
  snapbirth=getSub(trackId)['SnapshotIndexOfBirth']
  for isnap in range(snapbirth, MaxSnap+1):
	  f=h5py.File(rootdir+'/SubSnap_%03d.hdf5'%isnap,'r')
	  sat=f['Subhalos'][trackId]
	  HostTrackId=f['Membership/GroupedTrackIds'][sat['HostHaloId']][0]
	  hosttrackId.append(HostTrackId)
	  f.close()
  return np.array(hosttrackId)

f=h5py.File(rootdir+'SubSnap_%03d.hdf5'%MaxSnap,'r')
s=f['Subhalos'][...]
s.sort(order='TrackId')
selectionfunc=(s['SnapshotIndexOfDeath']-s['SnapshotIndexOfLastMaxMass']>50)&(s['HostHaloId']==0)
trackId=find(selectionfunc)[1]

sat,Host,host,snaps=getTrackPair(trackId)
#Host=getTrack(host[-1]['TrackId'])

def plotPair(XProxy='ComovingMostBoundPosition'):
  xsat=sat[XProxy]
  #xhost=host[XProxy]
  xHost=Host[XProxy]
  plt.figure()
  plt.plot(xsat[:,0], xsat[:,1],'ro-')
  plt.plot(xHost[:,0], xHost[:,1],'gx-')
  plt.figure()
  plt.plot(xsat[:,0], xsat[:,2],'ro-')
  plt.plot(xHost[:,0], xHost[:,2],'gx-')

def distance(x,y, axis=1):
  return np.sqrt(np.sum((x-y)**2, axis=axis))

def norm(x, axis=-1):
  return np.sqrt(np.sum(x**2, axis=axis))

def plotOrbit(XProxy='ComovingMostBoundPosition', VProxy='PhysicalAverageVelocity'):
  xsat=sat[XProxy]
  xHost=Host[XProxy]
  X=sat[XProxy]-Host[XProxy]
  V=sat[VProxy]-Host[VProxy]
  r=norm(X)
  pot,pots=NFWPot(sat, Host, snaps[:,1], snaps[:,2], XProxy)
  k=0.5*norm(V)**2
  E=pot+k
  J=norm(np.cross(X, V)**2)
  VmaxSat=sat['VmaxPhysical']
  VmaxHost=Host['VmaxPhysical']
  msat=sat['Nbound'];
  mhost=Host['Nbound'];
  rhost=Host['R200CritComoving'];
  snapend=sat[-1]['SnapshotIndexOfDeath']
  plt.figure()
  t=snaps[:,0]
  plt.subplot(221)
  plt.plot(t, r/rhost,'-o')
  plt.yscale('log')
  plt.plot([snapend, snapend], plt.ylim(),'k--')
  plt.ylabel('Separation')
  ax=plt.subplot(222)
  plt.plot(t, msat, 'r-')
  plt.plot(t, VmaxSat, 'r--')
  plt.plot(t, mhost, 'g-')
  plt.plot(t, VmaxHost, 'g--')
  plt.yscale('log')
  plt.plot([snapend, snapend], plt.ylim(),'k--')
  plt.ylabel('Mass and Vmax')
  ax.yaxis.set_label_position('right')
  plt.subplot(223)
  plt.plot(t, E/np.abs(E[-1]), '-')
  plt.plot([snapend, snapend], plt.ylim(),'k--')
  plt.ylabel('Energy')
  #plt.yscale('symlog', linthreshy=0.01)
  ax=plt.subplot(224)
  plt.plot(t, J, '-')
  plt.yscale('log')
  plt.plot([snapend, snapend], plt.ylim(),'k--')
  plt.ylabel('Angular momentum')
  ax.yaxis.set_label_position('right')
  
#def plotOrbit(XProxy='ComovingMostBoundPosition'):
  #nsnap=min(len(sat), len(Host))
  #xsat=sat[XProxy]
  #xHost=Host[XProxy]
  #r=distance(xsat[-nsnap:], xHost[-nsnap:])
  #msat=sat['Nbound'][-nsnap:];
  #mhost=Host['Nbound'][-nsnap:];
  #rhost=Host['R200CritComoving'][-nsnap:];
  #snaps=np.arange(nsnap)+MaxSnap-nsnap+1;
  #snapstart=MaxSnap-nsnap+1;
  #iend=sat[-1]['SnapshotIndexOfDeath']-snapstart
  ##iend=MaxSnap-snapstart;
  #plt.figure()
  #plt.plot(snaps[:iend], r[:iend]/rhost[:iend],'-o')
  #plt.figure()
  #plt.plot(snaps[:iend], msat[:iend], 'ro-')
  #plt.plot(snaps[:iend], mhost[:iend], 'gx-')
  #plt.yscale('log')

XProxy='ComovingAveragePosition';
#XProxy='ComovingMostBoundPosition';

#trackId=find(selectionfunc)[2]

#sat,host=getTrackPair(trackId)
#Host=getTrack(host[-1]['TrackId'])

#plotPair(XProxy)
plotOrbit('ComovingAveragePosition', 'PhysicalAverageVelocity')

plotOrbit('ComovingMostBoundPosition', 'PhysicalMostBoundVelocity')


import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from numpy.lib.recfunctions import append_fields

def PeriodicDistance(x,y, BoxSize, axis=-1):
  d=x-y
  d[d>BoxSize/2]=d[d>BoxSize/2]-BoxSize
  d[d<-BoxSize/2]=d[d<-BoxSize/2]+BoxSize
  return np.sqrt(np.sum(d**2, axis=axis))

def distance(x,y, axis=1):
  return np.sqrt(np.sum((x-y)**2, axis=axis))

def GetFileName(isnap, rootdir, ifile, singlefile=False, filetype='Sub'):
  if singlefile:
	return rootdir+'/'+filetype+'Snap_%03d.hdf5'%isnap
  else:
	return rootdir+'/'+filetype+'Snap_%03d.%d.hdf5'%(isnap,ifile)
  
def LoadSubhalos(isnap, rootdir):
  singlefile=False
  try:
	nfiles=h5py.File(GetFileName(isnap, rootdir, 0, False),'r')['NumberOfFiles'][...]
  except:
	singlefile=True
	nfiles=1
	
  subhalos=[]
  for i in xrange(nfiles):
	subfile=h5py.File(GetFileName(isnap, rootdir, i, singlefile), 'r')
	subhalos.append(subfile['Subhalos'][...])
	subfile.close()
  subhalos=np.hstack(subhalos)
  #subhalos.sort(order=['HostHaloId','Nbound'])
  return subhalos

def LoadParticles(isnap, rootdir, filetype='Sub'):
  singlefile=False
  try:
	nfiles=h5py.File(GetFileName(isnap, rootdir, 0, False, filetype),'r')['NumberOfFiles'][...]
  except:
	singlefile=True
	nfiles=1
	
  subhalos=[]
  for i in xrange(nfiles):
	subfile=h5py.File(GetFileName(isnap, rootdir, i, singlefile, filetype), 'r')
	subhalos.append(subfile[filetype+'haloParticles'][...])
	subfile.close()
  subhalos=np.hstack(subhalos)
  #subhalos.sort(order=['HostHaloId','Nbound'])
  return subhalos

def GetSub(trackId, isnap, rootdir):
  subhalos=LoadSubhalos(isnap, rootdir)
  return subhalos[subhalos['TrackId']==trackId]

def GetTrack(trackId, rootdir, MaxSnap):
  track=[];
  snaps=[]
  snapbirth=GetSub(trackId, MaxSnap, rootdir)['SnapshotIndexOfBirth']
  for isnap in range(snapbirth, MaxSnap+1):
	  s=GetSub(trackId, isnap, rootdir)
	  track.append(s)
	  snaps.append(isnap)
  return append_fields(np.array(track), 'Snapshot', np.array(snaps), usemask=False)

def GetScaleFactor(isnap, rootdir):
  try:
	a=h5py.File(GetFileName(isnap, rootdir, 0, False),'r')['ScaleFactor'][0]
  except:
	a=h5py.File(GetFileName(isnap, rootdir, 0, True),'r')['ScaleFactor'][0]

  return a
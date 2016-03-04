import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from numpy.lib.recfunctions import append_fields

def LoadSubhalos(isnap, rootdir):
  nfiles=h5py.File(rootdir+'SubSnap_%03d.%d.hdf5'%(isnap,0), 'r')['NumberOfFiles'][...]
  subhalos=[]
  for i in xrange(nfiles):
	subfile=h5py.File(rootdir+'SubSnap_%03d.%d.hdf5'%(isnap,i), 'r')
	subhalos.append(subfile['Subhalos'][...])
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
  snapbirth=getSub(trackId, MaxSnap, rootdir)['SnapshotIndexOfBirth']
  for isnap in range(snapbirth, MaxSnap+1):
	  s=getSub(trackId, isnap)
	  track.append(s)
	  snaps.append(isnap)
  return append_fields(np.array(track), 'Snapshot', np.array(snaps)).data
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from numpy.lib.recfunctions import append_fields

hbtdir="/gpfs/data/jvbq85/HBT/data/AqA5/subcat2/"
MaxSnap=64

def combineSubSnap(snapnum):
  subs=[]
  particles=[]
  NumberOfNewSubhalos=0
  NumberOfFakeHalos=0
  
  infile=h5py.File(hbtdir+'/SubSnap_%03d.0.hdf5'%snapnum,'r')
  outfile=h5py.File(hbtdir+'/SubSnap_%03d.hdf5'%snapnum, 'w')
  for item in ['SnapshotId',  'HubbleParam', 'ScaleFactor']: #'OmegaM0', 'OmegaLambda0',
	infile.copy(item, outfile['/'])
  NumberOfFiles=infile['NumberOfFiles'][...]
  infile.close()
  
  for ifile in xrange(NumberOfFiles):
	infile=h5py.File(hbtdir+'/SubSnap_%03d.%d.hdf5'%(snapnum, ifile), 'r')
	NumberOfNewSubhalos+=infile['NumberOfNewSubhalos'][0]
	NumberOfFakeHalos+=infile['NumberOfFakeHalos'][0]
	subs.append(infile['Subhalos'][...])
	particles.append(infile['SubhaloParticles'][...])

  subs=np.hstack(subs)
  order=np.argsort(subs, order='TrackId')
  outfile.create_dataset('Subhalos', data=subs[order])
  outfile.close()
  
  outfile=h5py.File(hbtdir+'/SubParticles_%03d.hdf5'%snapnum, 'w')
  particles=np.hstack(particles)
  dt = h5py.special_dtype(vlen=particles[0].dtype)
  dset=outfile.create_dataset('SubhaloParticles', particles.shape, dtype=dt)
  for i in xrange(len(particles)):
	dset[i]=particles[order[i]]
  outfile.close()

combineSubSnap(20)
#for isnap in xrange(MaxSnap):
  #combineSubSnap(isnap)
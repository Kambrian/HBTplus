import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from numpy.lib.recfunctions import append_fields

hbtdir="/gpfs/data/jvbq85/HBT/data/Millennium2/subcat2/"
mp=6.88e-4

def combineSubSnap(snapnum):
  subs=[]
  particles=[]
  nests=[]
  NumberOfNewSubhalos=0
  NumberOfFakeHalos=0
  
  infile=h5py.File(hbtdir+'/%03d/SubSnap_%03d.0.hdf5'%(snapnum,snapnum),'r')
  outfile=h5py.File(hbtdir+'/SubSnap_%03d.hdf5'%snapnum, 'w')
  infile.copy('Cosmology', outfile['/'])
  infile.copy('SnapshotId', outfile['/'])
  NumberOfFiles=infile['NumberOfFiles'][...]
  infile.close()
  
  outfile.create_dataset('Cosmology/ParticleMass', data=mp, shape=(1,), dtype='float32')
  for ifile in xrange(NumberOfFiles):
	infile=h5py.File(hbtdir+'/%03d/SubSnap_%03d.%d.hdf5'%(snapnum, snapnum, ifile), 'r')
	NumberOfNewSubhalos+=infile['NumberOfNewSubhalos'][0]
	NumberOfFakeHalos+=infile['NumberOfFakeHalos'][0]
	subs.append(infile['Subhalos'][...])
	particles.append(infile['SubhaloParticles'][...])
	nests.append(infile['NestedSubhalos'][...])
	infile.close()

  subs=np.hstack(subs)
  order=np.argsort(subs, order='TrackId')
  outfile.create_dataset('Subhalos', data=subs[order])
  outfile.create_dataset('NumberOfNewSubhalos', data=NumberOfNewSubhalos, shape=(1,))
  outfile.create_dataset('NumberOfFakeHalos', data=NumberOfFakeHalos, shape=(1,))
  
  write_vlen_data(particles, order, "SubhaloParticles", outfile)
  write_vlen_data(nests, order, "NestedSubhalos", outfile)
  outfile['NestedSubhalos'].attrs.create("Comment","List of the TrackIds of first-level sub-subhaloes within each subhalo.")
  
  outfile.close()
  
def write_vlen_data(data, order, name, outfile):
  print name, '....'
  sys.stdout.flush()
  data=np.hstack(data)[order]
  dt = h5py.special_dtype(vlen=data[0].dtype)
  dset=outfile.create_dataset(name, data.shape, dtype=dt, data=data)
  #for i in xrange(len(data)):
	#dset[i]=data[order[i]]

import shutil
#combineSubSnap(20)
for isnap in xrange(4, 68):
  print isnap
  sys.stdout.flush()
  combineSubSnap(isnap)
  #shutil.rmtree('%03d'%isnap)	
  
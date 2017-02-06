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
  nests=[]
  particleproperties=[]
  bindingenergy=[]
  
  NumberOfNewSubhalos=0
  NumberOfFakeHalos=0
  
  infile=h5py.File(hbtdir+'%03d/SubSnap_%03d.0.hdf5'%(snapnum,snapnum),'r')
  outfile=h5py.File(hbtdir+'/SubSnap_%03d.hdf5'%snapnum, 'w')
  infile.copy('Cosmology', outfile['/'])
  infile.copy('SnapshotId', outfile['/'])
  NumberOfFiles=infile['NumberOfFiles'][...]
  infile.close()
  
  for ifile in xrange(NumberOfFiles):
	infile=h5py.File(hbtdir+'%03d/SubSnap_%03d.%d.hdf5'%(snapnum, snapnum, ifile), 'r')
	NumberOfNewSubhalos+=infile['NumberOfNewSubhalos'][0]
	NumberOfFakeHalos+=infile['NumberOfFakeHalos'][0]
	subs.append(infile['Subhalos'][...])
	particles.append(infile['SubhaloParticles'][...])
	nests.append(infile['NestedSubhalos'][...])
	if 'ParticleProperties' in infile:
	  particleproperties.append(infile['ParticleProperties'][...])
	if 'BindingEnergies' in infile:
	  bindingenergy.append(infile['BindingEnergies'][...])

  outfile.create_dataset('NumberOfNewSubhalos', data=NumberOfNewSubhalos)
  outfile.create_dataset('NumberOfFakeHalos', data=NumberOfFakeHalos)
  subs=np.hstack(subs)
  order=np.argsort(subs, order='TrackId')
  outfile.create_dataset('Subhalos', data=subs[order])
  #outfile.close()
  #outfile=h5py.File(hbtdir+'/SubParticles_%03d.hdf5'%snapnum, 'w')
  
  write_vlen_data(particles, order, 'SubhaloParticles', outfile)
  write_vlen_data(nests, order, 'NestedSubhalos', outfile)
  if len(particleproperties):
	write_vlen_data(particleproperties, order, 'ParticleProperties', outfile)
  if len(bindingenergy):
	write_vlen_data(bindingenergy, order, 'BindingEnergies', outfile)
          
  outfile.close()
  
def write_vlen_data(data, order, name, outfile):
  data=np.hstack(data)
  dt = h5py.special_dtype(vlen=data[0].dtype)
  dset=outfile.create_dataset(name, data.shape, dtype=dt)
  for i in xrange(len(data)):
	dset[i]=data[order[i]]

combineSubSnap(20)
#for isnap in xrange(MaxSnap):
  #combineSubSnap(isnap)
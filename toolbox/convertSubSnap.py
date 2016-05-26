'''to convert subhalo files from (the more advanced) structured array to plain array datasets for each property separately.
will also convert the particle list into separate arrays of particles (or a combined single array of particles if combine_particlelist=True).
input: 
   inputdir: directory for existing subhalo files
   snapnum: snapshot number to convert
output:
   converted subhalo files under inputdir/converted/.
'''

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys,os
from numpy.lib.recfunctions import append_fields

inputdir=sys.argv[1]
snapnum=int(sys.argv[2])
#inputdir="/gpfs/data/jvbq85/HBT/data/AqA5/subcat2/"

outputdir=inputdir+'/converted/'
try:
  os.mkdir(outputdir)
except:
  pass

combine_particlelist=False #set to true to combine the particles into a single array

def convertSubSnap(snapnum):
  subs=[]
  particles=[]
  
  #copy header
  infile=h5py.File(inputdir+'/SubSnap_%03d.0.hdf5'%snapnum,'r')
  outfile=h5py.File(outputdir+'/SubSnap_%03d.hdf5'%snapnum, 'w')
  for item in ['SnapshotId',  'HubbleParam', 'ScaleFactor']: #'OmegaM0', 'OmegaLambda0',
	infile.copy(item, outfile['/'])
  NumberOfFiles=infile['NumberOfFiles'][...]
  infile.close()
  
  #load
  for ifile in xrange(NumberOfFiles):
	infile=h5py.File(inputdir+'/SubSnap_%03d.%d.hdf5'%(snapnum, ifile), 'r')
	subs.append(infile['Subhalos'][...])
	particles.append(infile['SubhaloParticles'][...])

  #save subhalo properties
  subs=np.hstack(subs)
  order=np.argsort(subs, order='TrackId')
  subs=subs[order];
  for field in subs.dtype.names:
	outfile.create_dataset('/Subhalos/'+field, data=subs[field])
  outfile.close()
  
  #save particles
  outfile=h5py.File(outputdir+'/SubParticles_%03d.hdf5'%snapnum, 'w')
  particles=np.hstack(particles)[order]
  if combine_particlelist:
	particles=np.hstack(list(particles))
	outfile.create_dataset('SubhaloParticles', data=particles)
  else:
	for i,p in enumerate(particles):
	  outfile.create_dataset('Subhalo'+str(i), data=p)
  outfile.close()

convertSubSnap(snapnum)
#for isnap in xrange(MaxSnap):
  #convertSubSnap(isnap)
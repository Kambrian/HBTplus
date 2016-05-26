import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os.path, glob
from numpy.lib.recfunctions import append_fields

def PeriodicDistance(x,y, BoxSize, axis=-1):
  d=x-y
  d[d>BoxSize/2]=d[d>BoxSize/2]-BoxSize
  d[d<-BoxSize/2]=d[d<-BoxSize/2]+BoxSize
  return np.sqrt(np.sum(d**2, axis=axis))

def distance(x,y, axis=1):
  return np.sqrt(np.sum((x-y)**2, axis=axis))

class ConfigReader:
  ''' class to read the config files '''
  def __init__(self, config_file):
	self.Options={}
	with open(config_file, 'r') as f:
	  for line in f:
		pair=line.lstrip().split("#",1)[0].split("[",1)[0].split()
		if len(pair)==2:
		  self.Options[pair[0]]=pair[1]
		elif len(pair)>2:
		  self.Options[pair[0]]=pair[1:]
  
  def __getitem__(self, index):
	return self.Options[index]
  
	
class HBTReader:
  ''' class to read HBT2 catalogue '''
  
  def __init__(self, config_file):
	''' initialize HBTReader according to parameters in the HBT configuration file, which is the file with which HBT was called. The config file can also be found inside the subhalo directory as VER***.param '''
	self.Options=ConfigReader(config_file).Options
	self.rootdir=self.Options['SubhaloPath']
	self.MaxSnap=int(self.Options['MaxSnapshotIndex'])
	self.BoxSize=float(self.Options['BoxSize'])
	self.Softening=float(self.Options['SofteningHalo'])

	lastfile=sorted(glob.glob(self.rootdir+'/'+'SubSnap_*.hdf5'))[-1]
	extension=lastfile.rsplit('SubSnap_')[1].split('.')
	MaxSnap=int(extension[0])
	if MaxSnap!=self.MaxSnap:
	  print "HBT run not finished yet, maxsnap %d found (expecting %d)"%(MaxSnap, self.MaxSnap)
	  self.MaxSnap=MaxSnap

	self.nfiles=0
	if len(extension)==3:
	  self.nfiles=int(extension[1])+1
	  print self.nfiles, "subfiles per snapshot"

	if 'MinSnapshotIndex' in self.Options:
	  self.MinSnap=int(self.Options['MinSnapshotIndex'])
	else:
	  self.MinSnap=0
	    
  def GetFileName(self, isnap, ifile=0, filetype='Sub'):
	if isnap<0:
	  isnap=self.MaxSnap+1+isnap
	if self.nfiles:
	  return self.rootdir+'/'+filetype+'Snap_%03d.%d.hdf5'%(isnap, ifile)
	else:
	  return self.rootdir+'/'+filetype+'Snap_%03d.hdf5'%(isnap)
  
  def LoadSubhalos(self, isnap=-1, fields=None, subindex=None):
	'''load subhalos from snapshot isnap (default =-1, means final snapshot; isnap<0 will count backward from final snapshot)
	if fields are given, only load the specified fields. otherwise load all fields.
	if subindex is given, only load the subhalo with that index. all fields of it will be loaded.
	
	...Note: subindex specifies the order of the subhalo in the file at the current snapshot, i.e., subhalo=AllSubhalo[subindex]. subindex==trackId for single file output, but subindex!=trackId for mpi multiple-file outputs. 
	'''
	subhalos=[]
	offset=0
	for i in xrange(max(self.nfiles,1)):
	  with h5py.File(self.GetFileName(isnap, i), 'r') as subfile:
		if subindex is None:
		  if fields is None:
			subhalos.append(subfile['Subhalos'][...])
		  else:
			subhalos.append(subfile['Subhalos'][fields])
		else:
		  nsub=subfile['Subhalos'].shape[0]
		  if offset+nsub>subindex:
			subhalos.append(subfile['Subhalos'][subindex-offset])
			break
		  offset+=nsub
	subhalos=np.hstack(subhalos)
	#subhalos.sort(order=['HostHaloId','Nbound'])
	return subhalos

  def LoadParticles(self, isnap=-1, subindex=None, filetype='Sub'):	  
	''' load subhalo particle list at snapshot isnap. 
	
	if subindex is given, only load subhalo of the given index (the order it appears in the file, subindex==trackId for single file output, but not for mpi multiple-file outputs). otherwise load all the subhaloes.
	
	default filetype='Sub' will load subhalo particles. set filetype='Src' to load source subhalo particles instead (for debugging purpose only).'''
	
	subhalos=[]
	offset=0
	for i in xrange(max(self.nfiles,1)):
	  with h5py.File(self.GetFileName(isnap,  i, filetype), 'r') as subfile:
		if subindex is None:
		  subhalos.append(subfile[filetype+'haloParticles'][...])
		else:
		  nsub=subfile['Subhalos'].shape[0]
		  if offset+nsub>subindex:
			subhalos.append(subfile[filetype+'haloParticles'][subindex-offset])
			break
		  offset+=nsub
	subhalos=np.hstack(subhalos)
	return subhalos

  def GetParticleProperties(self, subindex, isnap=-1):	  
	'''load subhalo particle properties for subhalo with index subindex (the order it appears in the file, subindex==trackId for single file output, but not for mpi multiple-file outputs)'''
	for i in xrange(max(self.nfiles,1)):
	  with h5py.File(self.GetFileName(isnap,  i, filetype), 'r') as subfile:
		nsub=subfile['Subhalos'].shape[0]
		if offset+nsub>subindex:
		  return subfile['ParticleProperties/Subhalo%d'%(subindex-offset)][...]
		offset+=nsub
	raise RuntimeError("subhalo %d not found"%subindex)
  
  def GetSub(self, trackId, isnap=-1):
	''' load a subhalo with the given trackId at snapshot isnap'''
	#subhalos=LoadSubhalos(isnap, rootdir)
	#return subhalos[subhalos['TrackId']==trackId]
	if self.nfiles:
	  subid=find(self.LoadSubhalos(isnap, 'TrackId')==trackId)[0]
	else:
	  subid=trackId
	return self.LoadSubhalos(isnap, subindex=subid)

  def GetTrack(self, trackId):
	''' load an entire track of the given trackId '''
	track=[];
	snaps=[]
	snapbirth=self.GetSub(trackId)['SnapshotIndexOfBirth']
	for isnap in range(snapbirth, self.MaxSnap+1):
		s=self.GetSub(trackId, isnap)
		track.append(s)
		snaps.append(isnap)
	return append_fields(np.array(track), 'Snapshot', np.array(snaps), usemask=False)

  def GetScaleFactor(self, isnap):
	try:
	  return h5py.File(self.GetFileName(isnap),'r')['Cosmology/ScaleFactor'][0]
	except:
	  return h5py.File(self.GetFileName(isnap),'r')['ScaleFactor'][0]


if __name__ == '__main__':
    import timeit
    #apostle=HBTReader('../configs/Apostle_S1_LR.conf')
    apostle=HBTReader('/cosma/home/jvbq85/data/HBT/data/apostle/S1_LR/subcat/VER1.8.1.param')
    #apostle=HBTReader('/cosma/home/jvbq85/data/HBT/data/MilliMill/subcat2_full/VER1.8.1.param')
    print(timeit.timeit("[apostle.LoadSubhalos(i, subindex=1) for i in range(10,apostle.MaxSnap)]", setup="from __main__ import apostle", number=1))
    print(timeit.timeit("[apostle.LoadSubhalos(i, 'Nbound') for i in range(10,apostle.MaxSnap)]", setup="from __main__ import apostle", number=1))
    print(timeit.timeit("apostle.LoadSubhalos(-1, 'Nbound')", setup="from __main__ import apostle", number=100))
    print(timeit.timeit("[apostle.LoadSubhalos(i) for i in range(10,apostle.MaxSnap)]", setup="from __main__ import apostle", number=1))
    print(timeit.timeit("apostle.GetTrack(12)", setup="from __main__ import apostle", number=1))
    print(timeit.timeit("apostle.GetTrack(103)", setup="from __main__ import apostle", number=1))
    
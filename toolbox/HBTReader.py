'''
class to read the HBT outputs.

To use it, initialize the reader with the parameter file under the subhalo directory, e.g.,

  from HBTReader import HBTReader

  reader=HBTReader('subcat')

  snapshotnumber=-1 # or 0~MaxSnap. -1 means last snapshot

  subs=reader.LoadSubhalos(snapshotnumber) #load all

  nbound=reader.LoadSubhalos(snapshotnumber, 'Nbound') #only Nbound

  sub2=reader.LoadSubhalos(snapshotnumber, subindex=2) #only subhalo 2

  track2=reader.GetTrack(2) #track 2

'''
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import h5py
import sys
import os.path, glob
from numpy.lib.recfunctions import append_fields
import numbers
from time import time

if sys.version_info[0] == 2:
    range = xrange


def PeriodicDistance(x,y, BoxSize, axis=-1):
    d=x-y
    d[d>BoxSize/2]=d[d>BoxSize/2]-BoxSize
    d[d<-BoxSize/2]=d[d<-BoxSize/2]+BoxSize
    return np.sqrt(np.sum(d**2, axis=axis))

def distance(x,y, axis=-1):
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

def get_hbt_snapnum(snapname):
    return int(snapname.rsplit('SubSnap_')[1].split('.')[0])


class HBTReader:
  ''' class to read HBT2 catalogue '''

  def __init__(self, subhalo_path):
    ''' initialize HBTReader to read data in subhalo_path. A parameter file must 
    exist there (Parameters.log dumped by HBT during runtime).'''
    config_file=subhalo_path+'/Parameters.log'
    self.Options=ConfigReader(config_file).Options
    self.rootdir=self.Options['SubhaloPath']
    self.MaxSnap=int(self.Options['MaxSnapshotIndex'])
    self.BoxSize=float(self.Options['BoxSize'])
    self.Softening=float(self.Options['SofteningHalo'])
    # due to recent changes in the cosma machines, things may have
    # moved around
    if not os.path.isdir(self.rootdir):
        self.rootdir = self.rootdir.replace(
            '/gpfs/data/jvbq85', '/cosma/home/durham/jvbq85/data')

    try:
      lastfile=sorted(glob.glob(self.rootdir+'/SubSnap_*.hdf5'),
                      key=get_hbt_snapnum)[-1]
    except:
      lastfile=sorted(glob.glob(self.rootdir+'/*/SubSnap_*.hdf5'),
                      key=get_hbt_snapnum)[-1]
    extension=lastfile.rsplit('SubSnap_')[1].split('.')
    MaxSnap=int(extension[0])
    if MaxSnap!=self.MaxSnap:
      print("HBT run not finished yet, maxsnap %d found (expecting %d)" \
            %(MaxSnap, self.MaxSnap))
      self.MaxSnap=MaxSnap

    self.nfiles=0
    if len(extension)==3:
      self.nfiles=len(glob.glob(
        self.rootdir+'/%03d'%MaxSnap+'/SubSnap_%03d.*.hdf5'%MaxSnap))
      print(self.nfiles, "subfiles per snapshot")

    if 'MinSnapshotIndex' in self.Options:
      self.MinSnap=int(self.Options['MinSnapshotIndex'])
    else:
      self.MinSnap=0

    try:
      with self.Open(-1) as f:
        self.ParticleMass=f['/Cosmology/ParticleMass'][0]
    except:
      #print("Info: fail to get ParticleMass.")
      pass

    with self.Open(-1) as f:
      self.OmegaM0=f['/Cosmology/OmegaM0'][0]
      self.OmegaLambda0=f['/Cosmology/OmegaLambda0'][0]

  def Snapshots(self):
    return np.arange(self.MinSnap, self.MaxSnap+1)

  def GetFileName(self, isnap, ifile=0, filetype='Sub'):
    if isnap<0:
      isnap=self.MaxSnap+1+isnap
    if self.nfiles:
      return self.rootdir+'/%03d/'%isnap+filetype+'Snap_%03d.%d.hdf5' \
            %(isnap, ifile)
    else:
      return self.rootdir+'/'+filetype+'Snap_%03d.hdf5'%(isnap)

  def Open(self, isnap, ifile=0, filetype='Sub', mode='r'):
    '''return the opened hdf5 file'''
    return h5py.File(self.GetFileName(isnap, ifile, filetype), mode)
  
  def LoadNestedSubhalos(self, isnap=-1, selection=None):
    '''load the list of nested subhalo indices for each subhalo'''
    nests=[]
    for i in range(max(self.nfiles,1)):
      with h5py.File(self.GetFileName(isnap, i), 'r') as subfile:
        nests.extend(subfile['NestedSubhalos'][...])
    return np.array(nests)

  def LoadSubhalos(self, isnap=-1, selection=None, show_progress=False):
    '''load subhalos from snapshot isnap (default =-1, means final snapshot; 
    isnap<0 will count backward from final snapshot)

    `selection` can be a single field, a list of the field names or a single 
    subhalo index. e.g., selection=('Rank', 'Nbound') will load only the Rank 
    and Nbound fields of subhaloes. selection=3 will only load subhalo with 
    subindex 3. Default will load all fields of all subhaloes.

    ...Note: subindex specifies the order of the subhalo in the file at the 
    current snapshot, i.e., subhalo=AllSubhalo[subindex].    subindex==trackId 
    for single file output, but subindex!=trackId for mpi multiple-file outputs. 

    You can also use numpy slice for selection, e.g., selection=np.s_[:10, 
    'Rank','HostHaloId'] will select the 'Rank' and 'HostHaloId' of the first 10 
    subhaloes. You can also specify multiple subhaloes by passing a list of 
    (ordered) subindex, e.g., selection=((1,2,3),). However, currently only a 
    single subhalo can be specified for multiple-file hbt data (not restricted 
    for single-file data).

    '''
    subhalos=[]
    offset=0
    trans_index=False
    if selection is None:
      selection=np.s_[:]
    else:
      trans_index=isinstance(selection, numbers.Integral)

    if type(selection) is list:
      selection=tuple(selection)

    for i in range(max(self.nfiles,1)):
      if show_progress:
        sys.stdout.write(".")
        sys.stdout.flush()
      file = self.GetFileName(isnap, i)
      if not os.path.isfile(file):
        continue
      with h5py.File(file, 'r') as subfile:
        nsub=subfile['Subhalos'].shape[0]
        if nsub==0:
          continue
        if trans_index:
          if offset+nsub>selection:
            subhalos.append(subfile['Subhalos'][selection-offset])
            break
          offset+=nsub
        else:
          subhalos.append(subfile['Subhalos'][selection])
    if len(subhalos):
      subhalos=np.hstack(subhalos)
    else:
      subhalos=np.array(subhalos)
    if show_progress:
      print()
    #subhalos.sort(order=['HostHaloId','Nbound'])
    return subhalos

  def GetNumberOfSubhalos(self, isnap=-1):
    with h5py.File(self.GetFileName(isnap, 0),'r') as f:
        if self.nfiles:
            return f['TotalNumberOfSubhalosInAllFiles'][...]
        else:
            return f['Subhalos'].shape[0]

  def LoadParticles(self, isnap=-1, subindex=None, filetype='Sub'):	  
    ''' load subhalo particle list at snapshot isnap. 

    if subindex is given, only load subhalo of the given index (the order it 
    appears in the file, subindex==trackId for single file output, but not for 
    mpi multiple-file outputs). otherwise load all the subhaloes.

    default filetype='Sub' will load subhalo particles. set filetype='Src' to 
    load source subhalo particles instead (for debugging purpose only).'''

    subhalos=[]
    offset=0
    for i in range(max(self.nfiles,1)):
      with h5py.File(self.GetFileName(isnap,  i, filetype), 'r') as subfile:
        if subindex is None:
          subhalos.append(subfile[filetype+'haloParticles'][...])
        else:
          nsub=subfile[filetype+'haloParticles'].shape[0]
          if offset+nsub>subindex:
            subhalos.append(subfile[filetype+'haloParticles'][subindex-offset])
            break
          offset+=nsub
    subhalos=np.hstack(subhalos)
    return subhalos

  def GetParticleProperties(self, subindex, isnap=-1):	  
    '''load subhalo particle properties for subhalo with index subindex (the 
    order it appears in the file, subindex==trackId for single file output, but 
    not for mpi multiple-file outputs)'''

    offset=0
    for i in range(max(self.nfiles,1)):
      with h5py.File(self.GetFileName(isnap,  i), 'r') as subfile:
        nsub=subfile['Subhalos'].shape[0]
        if offset+nsub>subindex:
          #for compatibility with old data
          try:
            return subfile['ParticleProperties/Sub%d'%(subindex-offset)][...]
          except: 
            return subfile['ParticleProperties'][subindex-offset]
        offset+=nsub
    raise RuntimeError("subhalo %d not found"%subindex)
  
  def GetSub(self, trackId, isnap=-1):
    ''' load a subhalo with the given trackId at snapshot isnap'''
    #subhalos=LoadSubhalos(isnap, rootdir)
    #return subhalos[subhalos['TrackId']==trackId]
    if self.nfiles:
      #subid=find(self.LoadSubhalos(isnap, 'TrackId')==trackId)[0]
        TrackIds = self.LoadSubhalos(isnap, 'TrackId')
        subid = TrackIds[np.where(TrackIds == trackId)][0]
        #subid = where(self.LoadSubhalos(isnap, 'TrackId')==trackId)[0]
    else:
        subid = trackId
    return self.LoadSubhalos(isnap, subid)

  def GetTrackSnapshot(self, trackId, isnap, fields=None):
    """Get track information for a single snapshot"""
    s = self.GetSub(trackId, isnap)
    if fields is not None:
        return s[fields]
    return s

  def GetTrack(self, trackId, fields=None):
    ''' load an entire track of the given trackId '''
    snapbirth=self.GetSub(trackId)['SnapshotIndexOfBirth']
    if hasattr(snapbirth, '__iter__'):
        snapbirth = snapbirth[0]
    snaps = np.arange(snapbirth, self.MaxSnap+1, dtype=int)
    to = time()
    track = []
    scales = -1 * np.ones(snaps.size)
    for i, isnap in enumerate(snaps):
        if os.path.isfile(self.GetFileName(isnap)):
            track.append(self.GetTrackSnapshot(trackId, isnap, fields=fields))
            scales[i] = self.GetScaleFactor(isnap)
    track = np.array(track)
    return append_fields(
        track, ['Snapshot', 'ScaleFactor'], [snaps,scales], usemask=False)

  def GetScaleFactor(self, isnap):
    try:
      return h5py.File(self.GetFileName(isnap),'r')['Cosmology/ScaleFactor'][0]
    except:
      return h5py.File(self.GetFileName(isnap),'r')['ScaleFactor'][0]

  def GetScaleFactorDict(self):
    ''' return a dictionary that maps snapshot_index to ScaleFactor'''
    return dict([(i, self.GetScaleFactor(i))
                 for i in range(self.MinSnap, self.MaxSnap+1)])

  def GetExclusiveParticles(self, isnap=-1):
    '''return an exclusive set of particles for subhaloes at isnap, by assigning 
    duplicate particles to the lowest mass subhaloes'''
    OriginPart=self.LoadParticles(isnap)
    OriginPart=zip(range(len(OriginPart)),OriginPart)
    comp_mass=lambda x: len(x[1])
    OriginPart.sort(key=comp_mass)
    repo=set()
    NewPart=[]
    for i,p in OriginPart:
      if len(p)>1:
        p=set(p)
        p.difference_update(repo)
        repo.update(p)
      NewPart.append((i,list(p)))
    comp_id=lambda x: x[0]
    NewPart.sort(key=comp_id)
    NewPart=[x[1] for x in NewPart]
    return NewPart

if __name__ == '__main__':
    import timeit
    #apostle=HBTReader('../configs/Apostle_S1_LR.conf')
    apostle=HBTReader('/cosma/home/jvbq85/data/HBT/data/apostle/S1_LR/subcat/VER1.8.1.param')
    #apostle=HBTReader('/cosma/home/jvbq85/data/HBT/data/MilliMill/subcat2_full/VER1.8.1.param')
    print(timeit.timeit("[apostle.LoadSubhalos(i, 1) for i in range(10,apostle.MaxSnap)]", setup="from __main__ import apostle", number=1))
    #print(timeit.timeit("[apostle.LoadSubhalos(i, np.s_['Nbound','Rank']) for i in range(10,apostle.MaxSnap)]", setup="from __main__ import apostle,np", number=1))
    print(timeit.timeit("[apostle.LoadSubhalos(i, 'Nbound') for i in range(10,apostle.MaxSnap)]", setup="from __main__ import apostle", number=1))
    print(timeit.timeit("apostle.LoadSubhalos(-1, ('Nbound','Rank'))", setup="from __main__ import apostle", number=100))
    print(timeit.timeit("[apostle.LoadSubhalos(i) for i in range(10,apostle.MaxSnap)]", setup="from __main__ import apostle", number=1))
    print(timeit.timeit("apostle.GetTrack(12)", setup="from __main__ import apostle", number=1))
    print(timeit.timeit("apostle.GetTrack(103)", setup="from __main__ import apostle", number=1))
    

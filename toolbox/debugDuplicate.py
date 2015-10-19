import numpy as np
import matplotlib.pyplot as plt
import h5py

rootdir='/gpfs/data/jvbq85/HBT/data/AqA5/subcat2/'

def get_sub(snapid, subid=22):
  f=h5py.File(rootdir+'SubSnap_%03d.hdf5'%snapid,'r');
  data=f['/Subhalos'][subid]
  f.close()
  return data

def test_unique(snapid, subid=22):
  f=h5py.File(rootdir+'SrcSnap_%03d.hdf5'%snapid,'r');
  data=f['/SrchaloParticles'][subid][...]
  f.close()
  n1=len(data)
  n2=len(set(data))
  print n1,n2
  return data

test_unique(32)
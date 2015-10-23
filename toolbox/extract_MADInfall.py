import numpy as np
import matplotlib.pyplot as plt
import h5py

rootdir='/gpfs/data/jvbq85/HBT/data/MADHALOS/nfw.infall_full/subcat2/'
def getData(isnap):
  f=h5py.File(rootdir+'SubSnap_%03d.hdf5'%isnap,'r')
  s=f['Subhalos'][...]
  m=s[1]['Nbound']
  d=s[0]['ComovingPosition']-s[1]['ComovingPosition']
  d=np.sum(d**2)**0.5
  f.close()
  return d,m

data=[getData(i) for i in xrange(121)]
np.savetxt(rootdir+'dist-mass.txt', np.array(data))
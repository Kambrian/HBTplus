import numpy as np
import matplotlib.pyplot as plt
import h5py

rootdir='/gpfs/data/jvbq85/HBT/data/majormerger/dynamic/subcat2/'
def getData(isnap):
  f=h5py.File(rootdir+'SubSnap_%03d.hdf5'%isnap,'r')
  s=f['Subhalos'][...]
  s.sort(order='TrackId')
  m=s['Nbound']
  x=s['ComovingPosition']
  d=s[0]['ComovingPosition']-s[1]['ComovingPosition']
  d=np.sum(d**2)**0.5
  f.close()
  out=[]
  out.append(m[0])
  out.extend(x[0])
  out.append(m[1])
  out.extend(x[1])
  out.append(d)
  return out

data=np.array([getData(i) for i in xrange(101)])
np.savetxt(rootdir+'history.txt', data)

plt.plot(data[:,0]/1e6)
plt.plot(data[:,4]/1e6)
plt.figure();
plt.plot(data[:,-1])
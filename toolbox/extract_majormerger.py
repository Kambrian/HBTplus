import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
icat=sys.argv[1]
isnap=int(sys.argv[2])

rootdir='/gpfs/data/jvbq85/HBT/data/majormerger/dynamic/subcat%s/'%icat
def getData(isnap):
  f=h5py.File(rootdir+'SubSnap_%03d.0.hdf5'%isnap,'r')
  s=f['Subhalos'][...]
  f.close()
  f=h5py.File(rootdir+'SubSnap_%03d.1.hdf5'%isnap,'r')
  s1=f['Subhalos'][...]
  f.close()
  s=np.hstack([s,s1])
  s.sort(order='TrackId')
  m=s['Nbound']
  vmax=s['VmaxPhysical']
  x=s['ComovingMostBoundPosition']
  d=s[0]['ComovingMostBoundPosition']-s[1]['ComovingMostBoundPosition']
  d=np.sum(d**2)**0.5
  out=[]
  out.append(m[0])
  out.extend(x[0])
  out.append(m[1])
  out.extend(x[1])
  out.append(d)
  out.append(vmax[0])
  out.append(vmax[1])
  return out

data=np.array([getData(i) for i in xrange(isnap)])
np.savetxt(rootdir+'history.txt', data)
plt.figure();
plt.plot(data[:,0]/1e6)
plt.plot(data[:,4]/1e6)
plt.savefig(rootdir+'/figure_1.png')
plt.figure();
plt.plot(data[:,8])
plt.savefig(rootdir+'/figure_2.png')
plt.figure()
plt.plot(data[:,9])
plt.plot(data[:,10])
plt.savefig(rootdir+'/figure_3.png')
plt.show()
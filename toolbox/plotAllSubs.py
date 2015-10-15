import numpy as np
import matplotlib.pyplot as plt
import h5py

rootdir='/gpfs/data/jvbq85/HBT/data/AqA5/'
f=h5py.File(rootdir+'subcat2/SubSnap_127.hdf5','r')
s2=f['SubHalos'][...]
f.close()

s=np.loadtxt(rootdir+'subcat/anal/halo_127',skiprows=1)

t=np.linspace(0,np.pi*2.+0.1, 20);
def plot_circle(x,y,r,**kwargs):
  return plt.plot(r*np.cos(t)+x,r*np.sin(t)+y,**kwargs)

def plot_seq(x,y,s,**kwargs):
  for i in range(len(x)):
	h=plot_circle(x[i],y[i],s[i],**kwargs)
  return h

massTosize=lambda m: m**0.33/800

plt.figure()
flt=s[:,3]>100
h1,=plot_seq(s[flt,4], s[flt,5], massTosize(s[flt,3]), color='g', linestyle='--', alpha=0.8)

flt=s2['Nbound']>100	
h2,=plot_seq(s2['ComovingPosition'][flt,0],s2['ComovingPosition'][flt,1], massTosize(s2['Nbound'][flt]), color='r', alpha=0.5)

plt.axis('equal')
plt.legend((h1,h2),('HBT','HBT2'))
plt.xlabel('x')
plt.ylabel('y')
plt.axis([56.8,57.3,52.4,52.8])
plt.savefig(rootdir+'subcat2/image3.pdf')

plt.figure()
flt=s[:,3]>0
plt.hist(np.log10(s[flt,3]),50, histtype='step', color='g', alpha=0.8)
flt=s2['Nbound']>0
plt.hist(np.log10(s2['Nbound'][flt]),50, histtype='step', color='r', alpha=0.5)
plt.yscale('log')
plt.legend(('HBT','HBT2'))
plt.xlabel(r'$\log10($Mass/Particles$)$')
plt.ylabel('Count')
plt.savefig(rootdir+'subcat2/hist.pdf')
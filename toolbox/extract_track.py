import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys


basedir='/gpfs/data/jvbq85/HBT/data/AqA5/'
snapbegin=0
snapend=127
def getTrack(trackId,catdir='subcat3'):
  rootdir=basedir+catdir
  track=[];
  snaps=[]
  for isnap in range(snapbegin, snapend+1):
	try:
	  s=h5py.File(rootdir+'/SubSnap_%03d.hdf5'%isnap,'r')['Subhalos'][trackId]
	  track.append(s)
	  snaps.append(isnap)
	except:
	  pass
	
  return np.array(track), np.array(snaps)

maintrack_v2=getTrack(27, 'subcat2_OMP')#from 21
maintrack_v3=getTrack(15, 'subcat3')
maintrack_v4=getTrack(16, 'subcat4')
sattrack_v4=getTrack(1477,'subcat4')
sat2track_v4=getTrack(212, 'subcat4')
#sat2track_v4=getTrack(5,'subcat4')

sattrack_v2=getTrack(1268, 'subcat2_OMP')#from 36
sattrack_v3=getTrack(1469, 'subcat3')#from 37
#d1=np.sqrt(np.sum((main[27-127-1:-1, 2:5]-sat[27-127-1:-1, 2:5])**2, axis=1))
d1_v2=np.sqrt(np.sum((maintrack_v2[0]['ComovingPosition'][27-127-1:-1]-sat2track_v2[0]['ComovingPosition'][27-127-1:-1])**2, axis=1))
d1_v3=np.sqrt(np.sum((maintrack_v3[0]['ComovingPosition'][27-127-1:-1]-sat2track_v3[0]['ComovingPosition'][27-127-1:-1])**2, axis=1))
d1_v4=np.sqrt(np.sum((maintrack_v4[0]['ComovingPosition'][27-127-1:-1]-sat2track_v4[0]['ComovingPosition'][27-127-1:-1])**2, axis=1))
plt.figure()
#plt.plot(d1, 'k-')
plt.plot(d1_v2, 'r-.')
plt.plot(d1_v3, 'g--')
plt.plot(d1_v4, 'b:')


#sat2track_v2=getTrack(4, 'subcat2_OMP') #from 24
#sat2track_v3=getTrack(3, 'subcat3') #from 25
sat2track_v2=getTrack(57, 'subcat2_OMP') #from 24
sat2track_v3=getTrack(79, 'subcat3') #from 25
d2_v2=np.sqrt(np.sum((maintrack_v2[0]['ComovingPosition'][21-127-1:]-sat2track_v2[0]['ComovingPosition'][21-127-1:])**2, axis=1))
d2_v3=np.sqrt(np.sum((maintrack_v3[0]['ComovingPosition'][21-127-1:]-sat2track_v3[0]['ComovingPosition'][21-127-1:])**2, axis=1))
plt.figure()
plt.plot(d2_v2, 'ro')
plt.plot(d2_v3, 'g-')

plt.figure()
plt.plot(maintrack_v2[1], maintrack_v2[0]['Nbound'], 'r--')
plt.plot(maintrack_v3[1], maintrack_v3[0]['Nbound'], 'g-')
plt.plot(sattrack_v2[1], sattrack_v2[0]['Nbound'], 'r--')
plt.plot(sattrack_v3[1], sattrack_v3[0]['Nbound'], 'g-')

dims=[0,2]
x=sattrack_v2[0]['ComovingPosition']
plt.plot(x[:,dims[0]], x[:,dims[1]], 'r-')
x=sattrack_v3[0]['ComovingPosition']
plt.plot(x[:,dims[0]], x[:,dims[1]], 'go')

x


from matplotlib import animation

track=[]
colors='rgb'

def get_data(i, track):
    if i>=track[1][0]:
	  t=track[0][i-track[1][0]]
	  return [t['ComovingPosition'][0],t['ComovingPosition'][1],t['Nbound']**0.33*10]
    else:
	  return []
    
def get_frames():
  frames=[]
  for snap in xrange(127):
    data=[]
    for i in xrange(len(tracks)):
	  d=get_data(snap, tracks[i])
	  if len(d)!=0:
		d.extend(colors[i])
		data.append(tuple(d))
    if len(data)>0:
	  data=np.array(data,dtype=[('x','f4'),('y','f4'),('s','f4'),('c','S1')])
	  #plt.scatter(data[0], data[1])
	  sca=plt.scatter(data['x'], data['y'], s=data['s'], c=data['c'])
	  sca.set_animated(1)
	  frames.append([sca])
  return frames

tracks=[maintrack_v2, sattrack_v2, sat2track_v2]
fig = plt.figure()
frames=get_frames()
ani = animation.ArtistAnimation(fig, frames, interval=20, blit=True,repeat_delay=1000)  

tracks=[maintrack_v3, sattrack_v3, sat2track_v3]
fig = plt.figure()
frames2=get_frames()
ani = animation.ArtistAnimation(fig, frames2, interval=20, blit=True,repeat_delay=1000)  

''' print snapshot list for apostle/EAGLE data

input snapshot path and the type to list: snap, snip or all (=snap+snip)
'''
import sys,os,glob

rootdir=sys.argv[1]
listtype=sys.argv[2] # snap, snip, all
#rootdir='/cosma/home/jvbq85/data/HBT/data/apostle/S1_HR/'

os.chdir(rootdir)
snaplist=sorted(glob.glob('snapshot_*'))
sniplist=sorted(glob.glob('snipshot_*'))

def redshift(snapname):
  return float(snapname.split('_z')[1].replace('p','.'))
  
comlist=sorted(snaplist+sniplist, key=redshift, reverse=True)
outlist={'snap':snaplist, 'snip':sniplist, 'all': comlist}

for x in outlist[listtype]:
  print x
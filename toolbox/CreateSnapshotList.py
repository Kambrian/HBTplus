''' print snapshot list for apostle/EAGLE data

input snapshot path and the type to list: snap, snip or all (=snap+snip)
'''
import sys,os,glob
import numpy as np

rootdir=sys.argv[1]
listtype=sys.argv[2] # snap, snip, all
#rootdir='/cosma/home/jvbq85/data/HBT/data/apostle/S1_HR/'

def redshift(snapname):
  return float(snapname.split('_z')[1].replace('p','.'))

def clear_dup(snaps):
  '''clear duplicate redshifts'''
  return list(np.array(snaps)[np.unique([redshift(s) for s in snaps], return_index=True)[1]])

def clear_empty(snaps):
  '''clear empty dirs'''
  return [x for x in snaps if os.listdir(x)!=[]]

os.chdir(rootdir)
snaplist=clear_empty(sorted(glob.glob('snapshot_*')))
n0=len(snaplist)
snaplist=clear_dup(snaplist)
sniplist=clear_empty(sorted(glob.glob('snipshot_*')))
n1=len(sniplist)
sniplist=clear_dup(sniplist)

comlist=sorted(list(snaplist)+list(sniplist), key=redshift, reverse=True)
n2=len(comlist)
comlist=clear_dup(comlist)
sys.stderr.write("%d duplicates cleared\n"%(n0+n1+n2-len(snaplist)-len(sniplist)-len(comlist)))

outlist={'snap':snaplist, 'snip':sniplist, 'all': comlist}

for x in outlist[listtype]:
  print x
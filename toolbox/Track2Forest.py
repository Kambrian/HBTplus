'''group tracks into forests of connected trees. 
useful for decomposing the dataset into indepedent domains for semi-analytical galaxy formation models.

usage: 

    python Track2Forest.py [subcat_path] <outfile>

If outputfile is not specified, then will write to subcat_path/TrackForest.hdf5


The output file contains three datasets: TrackForest, ForestSize, ForestOffset
The TrackIds and ForestIds can be accessed as TrackForest['TrackId'], TrackForest['ForestId']. 
The tracks in each forest can be found as  
        
        TrackForest['TrackId'][ForestOffset[forest_id]:ForestOffset[forest_id+1]]

'''

import numpy as np
import h5py
from HBTReader import HBTReader
from copy import deepcopy

def merge_sets(insets):
    '''merge intersecting sets. return list of merged sets'''
    outsets=deepcopy(insets) #copy to work on
    key2set=dict() #dict of known keys
    for setid,s in enumerate(outsets): #work on each set
        #collect known setids
        known_setids=set()
        for key in s:
            if key in key2set:
                known_setids.add(key2set[key])
        known_setids=list(known_setids)
        #merge known sets to s
        if len(known_setids)>0:
            known_setids.sort()
            for k in known_setids:
                s.update(outsets[k])
                outsets[k].clear()
        #update dict for keys in the merged set
        for key in s:
            key2set[key]=setid
    outsets=[s for s in outsets if s]
    return outsets

#datadir='/home/jxhan/data/simu/LCDM128/subcat'

def AssignForestIds(datadir):
    '''Assign a forestId to each trackId. Each forestId identifies a forest of all the tracks that are connected at some point.
    return TrackForest, a struct array with fields ['TrackId','ForestId']'''
    reader=HBTReader(datadir)

    Subhalos=reader.LoadSubhalos(-1, ['TrackId'])
    HostHistory=[set() for i in Subhalos] #list of hosttrackIds for each subhalo
    for isnap in xrange(reader.MinSnap, reader.MaxSnap+1, 1):
        Subhalos=reader.LoadSubhalos(isnap, ['TrackId','HostHaloId','Rank'])
        if Subhalos.size==0:
            continue
        #Subhalos.sort(order='TrackId')
        
        HostTrackId=reader.LoadHostTrackIds(isnap)
        HostTrackId[Subhalos['HostHaloId']<0]=Subhalos['TrackId'][Subhalos['HostHaloId']<0]
        #Subhalos=recfunctions.append_fields(Subhalos, 'HostTrackId', HostTrackId, usemask=False)
        
        for i,t in enumerate(Subhalos['TrackId']):
            HostHistory[t].add(HostTrackId[i])

    HostForests=merge_sets(HostHistory) #HostIds grouped into forests

    #build dictionary to lookup the forestId for each hostId
    Host2Forest=dict() 
    for forestid,forest in enumerate(HostForests):
        for hid in forest:
            Host2Forest[hid]=forestid
            
    #assign forestId to each track
    TrackForest=np.array([(t,Host2Forest[HostTrackId[i]]) for i,t in enumerate(Subhalos['TrackId'])], 
                        dtype=[('TrackId', Subhalos['TrackId'].dtype), ('ForestId', Subhalos['TrackId'].dtype)])
    
    return TrackForest

def SaveForest(TrackForest, outfile):
    
    TrackForest.sort(order='ForestId')
    
    ForestSize=np.bincount(TrackForest['ForestId'])
    ForestOffset=np.cumsum(ForestSize)
    ForestOffset=np.concatenate(([0],ForestOffset)) #The last ForestOffset is the total size of TrackForest.
    
    f=h5py.File(outfile, 'w')
    f.create_dataset('TrackForest', data=TrackForest)
    f.create_dataset('ForestSize', data=ForestSize)
    f.create_dataset('ForestOffset', data=ForestOffset)
    f.create_dataset('Comment', data="the tracks in each forest can be found as  #TrackForest['TrackId'][ForestOffset[forest_id]:ForestOffset[forest_id+1]]")
    f.close()


if __name__ == '__main__':
    import sys
    argc=len(sys.argv)
    if argc<2 or argc>3:
        print('Usage: ')
        print('     python '+sys.argv[0]+' [subcat_path] <outputfile> ')
        print('     Aggregate connected tracks into forests for subcat at subcat_path, and write the forest decomposition to outfile in hdf5 format')
        print('     If outputfile is not specified, then will write to subcat_path/TrackForest.hdf5\n')
        sys.exit()
    
    subcatdir=sys.argv[1]
    outfile=subcatdir+'/TrackForest.hdf5'
    if argc>2:
        outfile=sys.argv[2]


    TrackForest=AssignForestIds(subcatdir)
    SaveForest(TrackForest, outfile)

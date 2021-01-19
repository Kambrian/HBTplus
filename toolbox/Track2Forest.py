'''group tracks into forests of connected trees. 
useful for decomposing the dataset into indepedent domains for semi-analytical galaxy formation models.

usage: 

  from Track2Forest import AssignForestIds
  TrackForest=AssignForestIds(subcatdir)

#the TrackIds and ForestIds can now be accessed as TrackForest['TrackIds'], TrackForest['ForestIds']. You can also sort according to one of them, e.g.,
 
  TrackForest.sort(order='ForestId')

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
        
        try:
            with reader.Open(isnap) as f:
                GroupedTrackIds=f['Membership/GroupedTrackIds'][...]
            HostTrackId=np.array([GroupedTrackIds[s][0] for s in Subhalos['HostHaloId']])
        except: # alternative method when Membership is not available
            Host2Track=Subhalos[(Subhalos['Rank']==0)&(Subhalos['HostHaloId']>=0)][['HostHaloId','TrackId']]
            Host2Track=dict(zip(Host2Track['HostHaloId'], Host2Track['TrackId']))
            Host2Track.update({-1:-1})
            HostTrackId=np.array([Host2Track[s] for s in Subhalos['HostHaloId']])
        
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


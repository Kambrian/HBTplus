'''convert HBT tree to LHaloTree

The mergers are defined as follows in this version:
1) when SinkTrackId is present, the subhalo merges to SinkTrackId at SnapshotIndexOfDeath (which is no later than the sink snapshot. For orphans the detected sink could be after SnapshotIndexOfDeath.)
2) when SinkTrackId<0 and SnapshotIndexOfDeath>=0, the subhalo merges to the central subhalo at SnapshotIndexOfDeath.

Dead subhalos are cleaned (removed) from the Tree file.
'''

import numpy as np
import h5py
import sys
sys.path.append('/home/jxhan/HBTplus/toolbox')
from HBTReader import HBTReader
from copy import deepcopy
from collections import defaultdict
from numpy.lib import recfunctions

def LoadGroupMass(groupcatdir, nfiles_group, isnap):
    Group_M_Crit200=[]
    Group_M_Mean200=[]
    Group_M_TopHat=[]
    for ifile in range(nfiles_group):
        with h5py.File(groupcatdir+'/groups_%03d/groups_%03d.%i.hdf5'%(isnap, isnap, ifile), 'r') as f:
            Group_M_Crit200.append(f['Group/Group_M_Crit200'][...])
            Group_M_Mean200.append(f['Group/Group_M_Mean200'][...])
            Group_M_TopHat.append(f['Group/Group_M_TopHat200'][...])
    Group_M_Crit200=np.concatenate(Group_M_Crit200)
    Group_M_Mean200=np.concatenate(Group_M_Mean200)
    Group_M_TopHat=np.concatenate(Group_M_TopHat)
    masses=np.array([Group_M_Crit200,Group_M_Mean200,Group_M_TopHat]).T
    return masses

def LoadGroupMassHBT(subcatdir, isnap):
    with h5py.File(subcatdir+'/HaloSize/HaloSize_%03d.hdf5'%isnap, 'r') as f:
        Group_M_Crit200=f['HostHalos']['M200Mean']
        Group_M_Mean200=f['HostHalos']['M200Crit']
        Group_M_TopHat=f['HostHalos']['MVir']
    masses=np.array([Group_M_Crit200,Group_M_Mean200,Group_M_TopHat]).T
    return masses

class SnapAll_t:
    def __init__(self, subcatdir, groupcatdir, nfiles_group=0, groupsize_filetype='HBT'):
        '''load all the subhalo snapshots into Snapshots[] array
        groupcatdir tells the path to find FoF virial size files.
        groupsize_filetype can be HBT (if using HaloSize files generated by HBT) or Illustris (for illustris data)
        '''
        self.reader=HBTReader(subcatdir)
        #load all the subhalo snapshots with the needed fields
        self.Snapshots=[]
        self.MinSnap=self.reader.MinSnap
        self.MaxSnap=self.reader.MaxSnap
        self.ScaleFactors=[self.reader.GetScaleFactor(i) for i in range(self.MinSnap, self.MaxSnap+1)]
        field_list=['TrackId', 'HostHaloId', 'Rank', 'SnapshotIndexOfBirth', 'SnapshotIndexOfDeath',
                    'SinkTrackId', 'MostBoundParticleId', 'Nbound',  'Mbound', 'RHalfComoving', 'VmaxPhysical',
                    'ComovingMostBoundPosition', 'PhysicalAverageVelocity', 'SpecificAngularMomentum',
                    'SpecificSelfPotentialEnergy', 'SpecificSelfKineticEnergy']
        for isnap in range(self.reader.MinSnap, self.reader.MaxSnap+1, 1):
                Subhalos=self.reader.LoadSubhalos(isnap, field_list)
                if groupsize_filetype=='HBT':
                    HostMasses=LoadGroupMassHBT(subcatdir, isnap)
                else:
                    HostMasses=LoadGroupMass(groupcatdir, nfiles_group, isnap)
                SubHostMass=np.zeros([len(Subhalos),3], 'float')
                selcen=(Subhalos['Rank']==0)&(Subhalos['HostHaloId']>=0)
                SubHostMass[selcen]=HostMasses[Subhalos[selcen]['HostHaloId']]
                #Subhalos.sort(order='TrackId')
                HostTrackId=self.reader.LoadHostTrackIds(isnap)
                HostTrackId[Subhalos['HostHaloId']<0]=Subhalos['TrackId'][Subhalos['HostHaloId']<0] #FirstHaloInFoFGroup
                HostTrackId=np.array(HostTrackId)
                NextTrackInFoF=np.zeros_like(HostTrackId)-1
                with self.reader.OpenFile(isnap) as f:
                    TrackIdGroups = f['Membership/GroupedTrackIds'][...]
                for g in TrackIdGroups[:-1]:
                    if len(g)>1:
                        NextTrackInFoF[g[:-1]]=g[1:]
                Subhalos=recfunctions.append_fields(Subhalos, ['HostTrackId', 'NextTrackInFoF', 'HostM200Crit', 'HostM200Mean', 'HostMvir'],
                                                    [HostTrackId, NextTrackInFoF, SubHostMass.T[0], SubHostMass.T[1], SubHostMass.T[2]],
                                                    usemask=False)

                self.Snapshots.append(Subhalos)

    def GetTrack(self, trackId, fields=None):
        ''' load an entire track of the given trackId '''
        track = []
        snaps = []
        snapbirth = self.Snapshots[-1][trackId]['SnapshotIndexOfBirth']
        for isnap in range(snapbirth, self.MaxSnap+1):
            s = self.Snapshots[isnap-self.MinSnap][trackId]
            if fields is not None:
                s = s[fields]
            track.append(s)
            snaps.append(isnap)
        snaps=np.array(snaps)
        null_array=np.zeros_like(snaps)-1
        track=recfunctions.append_fields(np.array(track),
                                         ['SnapNum', 'NodeId', 'NextPro', 'FirstProTrack', 'FirstProSnap', 'DescendantTrack', 'DescendantSnap'],
                                         [snaps, null_array, null_array, null_array, null_array, null_array, null_array], usemask=False)
        return track

def GetDestiny(track):
    '''define the merger events for dead subhalos'''
    sub0=track[-1]
    SnapBirth=sub0['SnapshotIndexOfBirth']
    if sub0['SnapshotIndexOfDeath']>=0: #dead
        if sub0['SinkTrackId']>=0: #has sink
            DestTrackId=sub0['SinkTrackId']
            DestSnap=sub0['SnapshotIndexOfDeath'] #set to death or sink snap?
            #set to sink snap
            #for s in track[-1::-1]:
            #    if s['SinkTrackId']<0:
            #        break
            #    DestSnap=s['SnapNum']+1
        else: #no sink, only disrupted
            DestSnap=sub0['SnapshotIndexOfDeath'] #set to death snap
            DestTrackId=track[DestSnap-SnapBirth]['HostTrackId'] #set to CentralTrack at SnapDeath
            #DestTrackId=sub0['HostTrackId'] #alternatively, set to final CentralTrack
            #do we worry about those dead outside any normal halo?
            #if DestTrackId==sub0['TrackId']: #dead within itself
            #    DestTrackId=-1 #maybe this will be assigned automatically when filling the tree, so no need to do it.
    else: #alive
        DestTrackId=sub0['TrackId']
        DestSnap=sub0['SnapNum']

    return DestSnap,DestTrackId


class Tree_t:
    def __init__(self, forestId, TrackForest, ForestOffset, SnapDB):
        ''' collect tracks in the Forest specified by forestId, and merge them to a LTree'''
        self.TreeId=forestId
        self.TrackList=TrackForest['TrackId'][ForestOffset[forestId]:ForestOffset[forestId+1]]
        self.NumTracks=len(self.TrackList)
        self.Tracks, self.TrackMeta=self._fill(SnapDB)
        self._sort_tracks()
        self._clean()
        self.Tree=self._link()

    def _fill(self, SnapDB):
        Tracks=[SnapDB.GetTrack(tid) for tid in self.TrackList]
        TrackMeta=[]
        for track in Tracks:
            s,i=GetDestiny(track)
            TrackMeta.append((i,s,track[-1]['Nbound']))

        sub0=SnapDB.Snapshots[-1][0]
        TrackMeta=np.array(TrackMeta, dtype=[('DestTrackId', sub0['TrackId'].dtype),
                                    ('DestSnap', 'int32'), ('Nbound', sub0['Nbound'].dtype)])
        return Tracks, TrackMeta

    def _sort_tracks(self):
        #sort according to root, merger time and current size
        self.TrackMeta['DestSnap']*=-1 #to sort in descending order of DestSnap
        track_order=np.argsort(self.TrackMeta, order=('DestTrackId', 'DestSnap', 'Nbound'))
        self.TrackMeta['DestSnap']*=-1

        self.Tracks=[self.Tracks[i] for i in track_order]
        self.TrackMeta=self.TrackMeta[track_order]

    def _clean(self):
        #clean-up empty nodes
        for i in range(self.NumTracks):
            self.Tracks[i]=np.delete(self.Tracks[i], self.Tracks[i]['Nbound']<2)

    def _link(self):
        #fill in pro-desc info
        TrackDict=defaultdict(lambda :-1)
        nodeid=0
        for i in range(self.NumTracks):
            track=self.Tracks[i]
            tid=track[0]['TrackId']
            track['FirstProTrack']=track['TrackId']
            track['FirstProSnap'][1:]=track['SnapNum'][:-1]
            track['FirstProTrack'][0]=-1
            track['FirstProSnap'][0]=-1
            track['DescendantTrack']=track['TrackId']
            track['DescendantSnap'][:-1]=track['SnapNum'][1:]
            track['DescendantTrack'][-1]=self.TrackMeta[i]['DestTrackId']
            track['DescendantSnap'][-1]=self.TrackMeta[i]['DestSnap']
            if track['DescendantTrack'][-1]==track['TrackId'][-1]: #merge to itself, host track
                track['DescendantTrack'][-1]=-1
                track['DescendantSnap'][-1]=-1
            for j in range(len(track)):
                track[j]['NodeId']=nodeid
                TrackDict[(track[j]['TrackId'], track[j]['SnapNum'])]=nodeid #map (trackid, snapnum) to nodeid
                nodeid+=1

        TrackDict[(-1, -1)]=-1

        #join to tree
        Tree=np.hstack(self.Tracks)

        #update links
        for node in Tree:
            #reuse these columns to store node ids
            node['FirstProTrack']=TrackDict[(node['FirstProTrack'], node['FirstProSnap'])] #FirstProgenitor
            node['DescendantTrack']=TrackDict[(node['DescendantTrack'], node['DescendantSnap'])] #Descendant
            node['HostTrackId']=TrackDict[(node['HostTrackId'], node['SnapNum'])] #FirstHaloInFoFGroup
            node['NextTrackInFoF']=TrackDict[(node['NextTrackInFoF'], node['SnapNum'])] #NextHaloInFoFGroup

        #build next progenitor link
        for node in Tree:
            if node['DescendantTrack']>=0:
                proid=Tree[node['DescendantTrack']]['FirstProTrack']
                if proid!=node['NodeId']: #not main branch, need to attach as nextpro
                    nextpro=Tree[proid]['NextPro']
                    while nextpro>=0: #walk to end of pro chain
                        proid=nextpro
                        nextpro=Tree[proid]['NextPro']
                    Tree[proid]['NextPro']=node['NodeId'] #attach

        return Tree

    def save(self, filename):
        G=43007.1
        tree=self.Tree
        with h5py.File(filename, 'a') as f:
            g=f.create_group('Tree%d'%self.TreeId)
            g.create_dataset('SubhaloNumber', data=tree['TrackId'])
            g.create_dataset('Descendant', data=tree['DescendantTrack'])
            g.create_dataset('FirstProgenitor', data=tree['FirstProTrack'])
            g.create_dataset('NextProgenitor', data=tree['NextPro'])
            g.create_dataset('FirstHaloInFoFGroup', data=tree['HostTrackId'])
            g.create_dataset('NextHaloInFoFGroup', data=tree['NextTrackInFoF'])
            #g.create_dataset('FileNr', data=np.zeros(len(tree), 'int'))
            g.create_dataset('SnapNum', data=tree['SnapNum'])
            g.create_dataset('Group_M_Crit200', data=tree['HostM200Crit'])
            g.create_dataset('Group_M_Mean200', data=tree['HostM200Mean'])
            g.create_dataset('Group_M_TopHat200', data=tree['HostMvir'])
            g.create_dataset('SubhaloHalfMassRad', data=tree['RHalfComoving'])
            g.create_dataset('SubhaloIDMostBound', data=tree['MostBoundParticleId'])
            g.create_dataset('SubhaloLen', data=tree['Nbound'])
            #g.create_dataset('SubhaloLenType', data=tree['NboundType'])
            g.create_dataset('SubhaloMass', data=tree['Mbound'])
            #g.create_dataset('SubhaloMassType', data=tree['MboundType'])
            g.create_dataset('SubhaloPos', data=tree['ComovingMostBoundPosition'])
            g.create_dataset('SubhaloVel', data=tree['PhysicalAverageVelocity'])
            spin=np.sqrt((tree['SpecificAngularMomentum']**2).sum())*np.sqrt(
                np.abs(tree['SpecificSelfPotentialEnergy']+0.5*tree['SpecificSelfKineticEnergy'])
                )/G/tree['Mbound'] #peebles spin
            g.create_dataset('SubhaloSpin', data=spin)
            g.create_dataset('SubhaloVmax', data=tree['VmaxPhysical'])
            g.create_dataset('SubhaloVelDisp', data=np.sqrt(2.*tree['SpecificSelfKineticEnergy']))
            #g.close()

subcatdir='/home/cossim/IllustrisTNG/Illustris-3/subcat'
groupcatdir='/home/cossim/IllustrisTNG/Illustris-3/output'
nfiles_group=2
groupsize_filetype='Illustris'

#subcatdir='/home/jxhan/data/simu/8213/subcatNew'
#groupcatdir='/home/jxhan/data/simu/8213/fof'
#nfiles_group=0
#groupsize_filetype='HBT'

tree_filename=subcatdir+'/LTree.hdf5'

with h5py.File(subcatdir+'/TrackForest.hdf5','r') as f:
    TrackForest=f['TrackForest'][...]
    ForestOffset=f['ForestOffset'][...]
#ForestSize=np.diff(ForestOffset)

SnapDB=SnapAll_t(subcatdir, groupcatdir, nfiles_group, groupsize_filetype)

TreeNHalos=[]
for forestId in range(len(ForestOffset)-1):
    tree=Tree_t(forestId, TrackForest, ForestOffset, SnapDB)
    tree.save(tree_filename)
    TreeNHalos.append(len(tree.Tree))

TreeNHalos=np.array(TreeNHalos)
TotNSubhalos=np.array([np.sum(s['Nbound']>1) for s in SnapDB.Snapshots])
redshifts=1./np.array(SnapDB.ScaleFactors)-1

with h5py.File(tree_filename, 'a') as f:
    g=f.create_group('Header')
    g.create_dataset('Redshifts', data=redshifts)
    g.create_dataset('TotNSubhalos', data=TotNSubhalos)
    g.create_dataset('TreeNHalos', data=TreeNHalos)
    g.create_dataset('FirstSnapshotNr', data=SnapDB.MinSnap)
    g.create_dataset('LastSnapshotNr', data=SnapDB.MaxSnap)
    g.create_dataset('NTreesPerFile', data=len(TreeNHalos))
    #g.create_dataset('ParticleMass', data=SnapDB.reader.ParticleMass)

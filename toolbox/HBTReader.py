"""HBTReader module.

Python module containing ``HBTReader`` class, used for analysis of  HBT+ halo
finder output.

Example:

.. code::python

    from HBTReader import HBTReader
    reader=HBTReader('subcat')
    snapshotnumber=-1                                    # or 0~MaxSnap. -1 means last snapshot
    subs=reader.LoadSubhalos(snapshotnumber)             # load all
    nbound=reader.LoadSubhalos(snapshotnumber, 'Nbound') # only Nbound
    sub2=reader.LoadSubhalos(snapshotnumber, subindex=2) # only subhalo 2
    track2=reader.GetTrack(2)                            # track 2
"""
#!/usr/bin/env python
import sys
import glob
import numbers
import logging
import numpy as np
import h5py
from numpy.lib.recfunctions import append_fields


def PeriodicDistance(x, y, BoxSize, axis=-1):
    d = x - y
    d[d > BoxSize / 2] = d[d > BoxSize / 2] - BoxSize
    d[d < -BoxSize / 2] = d[d < -BoxSize / 2] + BoxSize
    return np.sqrt(np.sum(d**2, axis=axis))


def distance(x, y, axis=-1):
    return np.sqrt(np.sum((x - y)**2, axis=axis))


class ConfigReader:
    """Class to read the config files (i.e. ``Parameters.log``).
    """

    def __init__(self, config_file):
        self.Options = {}
        with open(config_file, 'r') as file:
            for line in file:
                pair = line.lstrip().split("#", 1)[0].split("[", 1)[0].split()
                if len(pair) == 2:
                    self.Options[pair[0]] = pair[1]
                elif len(pair) > 2:
                    self.Options[pair[0]] = pair[1:]

    def __getitem__(self, index):
        return self.Options[index]


def get_hbt_snapnum(snapname):
    """Extracts snapshot number from a filename.

    Arguments:
        snapname (str): name of snapshot file

    Returns:
        (int): snapshot number
    """
    return int(snapname.rsplit('SubSnap_')[1].split('.')[0])


class HBTReader:
    """Class to read the HBTPlus outputs.

    To use it, initialize the reader with the directory in which
    ``Parameters.log`` is stored - it is written by HBT during runtime.

    Arguments:
        subhalo_path (str): directory with config files
    """

    def __init__(self, subhalo_path):
        config_file = subhalo_path + '/Parameters.log'
        self.Options = ConfigReader(config_file).Options
        self.rootdir = self.Options['SubhaloPath']
        self.MaxSnap = int(self.Options['MaxSnapshotIndex'])
        self.BoxSize = float(self.Options['BoxSize'])
        self.Softening = float(self.Options['SofteningHalo'])

        try:
            lastfile = sorted(
                glob.glob(self.rootdir + '/SubSnap_*.hdf5'),
                key=get_hbt_snapnum)[-1]
        except:
            lastfile = sorted(
                glob.glob(self.rootdir + '/*/SubSnap_*.hdf5'),
                key=get_hbt_snapnum)[-1]

        extension = lastfile.rsplit('SubSnap_')[1].split('.')
        MaxSnap = int(extension[0])
        if MaxSnap != self.MaxSnap:
            logging.info(
                "HBT run not finished yet, maxsnap %d found (expecting %d)" %
                (MaxSnap, self.MaxSnap))
            self.MaxSnap = MaxSnap

        self.nfiles = 0
        if len(extension) == 3:
            self.nfiles = len(
                glob.glob(self.rootdir + '/%03d' % MaxSnap +
                          '/SubSnap_%03d.*.hdf5' % MaxSnap))
            logging.info("%d subfiles per snapshot" % self.nfiles)

        if 'MinSnapshotIndex' in self.Options:
            self.MinSnap = int(self.Options['MinSnapshotIndex'])
        else:
            self.MinSnap = 0

        try:
            with self.OpenFile(-1) as file:
                self.ParticleMass = file['/Cosmology/ParticleMass'][0]
        except:
            logging.error("Failed to get ParticleMass")

        with self.OpenFile(-1) as file:
            self.OmegaM0 = file['/Cosmology/OmegaM0'][0]
            self.OmegaLambda0 = file['/Cosmology/OmegaLambda0'][0]

    def Snapshots(self):
        return np.arange(self.MinSnap, self.MaxSnap + 1)

    def GetFileName(self, isnap, ifile=0, filetype='Sub'):
        """Returns filename of an HBT snapshot

        Arguments:
            isnap (int): snapshot of the file
            ifile (int): (default=0) index for sub-snapshots
            filetype (str): (default='Sub') 'Src', 'Sub' or 'HaloSize'

        Returns:
            (str): HBT snaphost filename
        """
        if isnap < 0:
            isnap = self.MaxSnap + 1 + isnap
        if self.nfiles:
            return self.rootdir + '/%03d/' % isnap + filetype + 'Snap_%03d.%d.hdf5' % (
                isnap, ifile)
        else:
            if filetype == 'HaloSize':
                return self.rootdir + '/HaloSize/' + 'HaloSize_%03d.hdf5' % (
                    isnap)
            else:
                return self.rootdir + '/' + filetype + 'Snap_%03d.hdf5' % (
                    isnap)

    def OpenFile(self, isnap, ifile=0, filetype='Sub', mode='r'):
        """Opens HDF5 file.

        Arguments:
            isnap (int): snapshot of the file
            ifile (int): (default=0) index for sub-snapshots
            filetype (str): (default='Sub') 'Src', 'Sub' or 'HaloSize'
            mode (chr): (default='r') file handle mode

        Returns:
            (File): HDF5 HBT file handle
        """
        return h5py.File(self.GetFileName(isnap, ifile, filetype), mode)

    def LoadNestedSubhalos(self, isnap=-1, selection=None):
        """Load the list of nested subhalo indices for each subhalo

        Arguments:
            isnap (int): (default = -1) snapshot number

        Returns:
            (numpy.ndarray): array of nested indices
        """
        nests = []
        for i in xrange(max(self.nfiles, 1)):
            with self.OpenFile(isnap, i) as subfile:
                nests.extend(subfile['NestedSubhalos'][...])
        return np.array(nests)

    def LoadSubhalos(self, isnap=-1, selection=None, show_progress=False):
        """Load all subhaloes from a snapshot.

        .. Note::

            ``selection=('Rank', 'Nbound')`` will load only the Rank and Nbound
            fields of subhaloes; ``selection=3`` will only load subhalo with
            subindex 3; default will load all fields of all subhaloes.  You can
            also use ``numpy`` slice for selection, e.g. ``selection=np.s_[:10,
            'Rank','HostHaloId']`` will select the ``Rank`` and ``HostHaloId``
            of the first 10 subhaloes. You can also specify multiple subhaloes
            by passing a list of (ordered) subindex, e.g.,
            ``selection=((1,2,3),)``.  However, currently only a single subhalo
            can be specified for multiple-file HBT data (not restricted for
            single-file data).

        .. Note::

            Subindex specifies the order of the subhalo in the file at the current
            snapshot, i.e., ``subhalo=AllSubhalo[subindex]``.  ``subindex == trackId``
            for single file output, but ``subindex != trackId`` for mpi multiple-file
            outputs.

        Arguments:
            isnap (int): (default = -1) snapshot
            selection (numpy.s\_): (default = None) can be a single field, a list of
                the field names or a single subhalo index
            show_progress (bool): (default = False)
        """
        subhalos = []
        offset = 0
        trans_index = False
        if selection is None:
            selection = np.s_[:]
        else:
            trans_index = isinstance(selection, numbers.Integral)

        if type(selection) is list:
            selection = tuple(selection)

        for i in xrange(max(self.nfiles, 1)):
            if show_progress:
                sys.stdout.write(".")
                sys.stdout.flush()
            with self.OpenFile(isnap, i) as subfile:
                nsub = subfile['Subhalos'].shape[0]
                if nsub == 0:
                    continue
                if trans_index:
                    if offset + nsub > selection:
                        subhalos.append(
                            subfile['Subhalos'][selection - offset])
                        break
                    offset += nsub
                else:
                    subhalos.append(subfile['Subhalos'][selection])

        if len(subhalos):
            subhalos = np.hstack(subhalos)
        else:
            subhalos = np.array(subhalos)
        if show_progress:
            print("")
        #subhalos.sort(order=['HostHaloId','Nbound'])

        return subhalos

    def GetNumberOfSubhalos(self, isnap=-1):
        """Retunrs number of subhaloes in a snapshot.

        Arguments:
            isnap (int): (default = -1) snapshot number
        """
        with self.OpenFile(isnap) as file:
            if self.nfiles:
                return file['TotalNumberOfSubhalosInAllFiles'][...]
            else:
                return file['Subhalos'].shape[0]

    def LoadParticles(self, isnap=-1, subindex=None, filetype='Sub'):
        """Loads subhalo particle list at snapshot

        If ``subindex`` is given, only load subhalo of the given index (the order it
        appears in the file, subindex==trackId for single file output, but not for
        mpi multiple-file outputs). Otherwise loads all the subhaloes.

        Default filetype (``Sub``) will load subhalo particles. Filetype ``Src``
        loads source subhalo particles instead (for debugging purpose only).

        Arguments:
            isnap (int): (default=-1) snapshot number
            subindex (int): (default=None) index of a subhalo
            filetype (str): (default='Sub') HBT file type
        """

        subhalos = []
        offset = 0
        for i in xrange(max(self.nfiles, 1)):
            with self.OpenFile(isnap, i, filetype) as subfile:
                if subindex is None:
                    subhalos.append(subfile[filetype + 'haloParticles'][...])
                else:
                    nsub = subfile[filetype + 'haloParticles'].shape[0]
                    if offset + nsub > subindex:
                        subhalos.append(
                            subfile[filetype + 'haloParticles'][subindex
                                                                - offset])
                        break
                    offset += nsub
        subhalos = np.hstack(subhalos)
        return subhalos

    def GetParticleProperties(self, subindex, isnap=-1):
        """Returns subhalo particle properties for subhalo with index subindex.

        Values are returned in the order they appear in the file,
        ``subindex==trackId`` for single file output (but not for mpi
        multiple-file outputs)
        """
        offset = 0
        for i in xrange(max(self.nfiles, 1)):
            with self.OpenFile(isnap, i) as subfile:
                nsub = subfile['Subhalos'].shape[0]
                if offset + nsub > subindex:
                    try:
                        return subfile['ParticleProperties/Sub%d' % (
                            subindex - offset
                        )][...]  #for compatibility with old data
                    except:
                        return subfile['ParticleProperties'][subindex - offset]
                offset += nsub
        raise RuntimeError("subhalo %d not found" % subindex)

    def GetSub(self, trackId, isnap=-1):
        """Loads a subhalo with the given ``trackId`` at snapshot ``isnap``.
        """
        #subhalos=LoadSubhalos(isnap, rootdir)
        #return subhalos[subhalos['TrackId']==trackId]
        if self.nfiles:
            subids = self.LoadSubhalos(isnap, 'TrackId')
            subid = subids[subids == trackId][0]
        else:
            subid = trackId
        return self.LoadSubhalos(isnap, subid)

    def GetTrack(self, trackId, fields=None， MaxSnap=None):
        """Loads an entire track of the given ``trackId``.
        """
        track = []
        snaps = []
        scales = []
        snapbirth = self.GetSub(trackId)['SnapshotIndexOfBirth']

        if MaxSnap is None:
            isnap_range = range(snapbirth, self.MaxSnap + 1)
        else:
            isnap_range = range(snapbirth, MaxSnap + 1)

        for isnap in isnap_range:
            s = self.GetSub(trackId, isnap)
            a = self.GetScaleFactor(isnap)
            if fields is not None:
                s = s[fields]
            track.append(s)
            snaps.append(isnap)
            scales.append(a)
        return append_fields(
            np.array(track), ['Snapshot', 'ScaleFactor'],
            [np.array(snaps), np.array(scales)],
            usemask=False)

    def GetScaleFactor(self, isnap):
        """Reads scale factor at a given snapshot.
        """
        try:
            return self.OpenFile(isnap)['Cosmology/ScaleFactor'][0]
        except:
            return self.OpenFile(isnap)['ScaleFactor'][0]

    def GetScaleFactorDict(self):
        """Returns a dictionary that maps ``snapshot_index`` to ``ScaleFactor``.
        """
        return dict([(i, self.GetScaleFactor(i))
                     for i in range(self.MinSnap, self.MaxSnap + 1)])

    def GetExclusiveParticles(self, isnap=-1):
        """Loads an exclusive set of particles for subhaloes at ``isnap``

        Duplicate particles are assigned to the lowest mass subhaloes.
        """
        OriginPart = self.LoadParticles(isnap)
        OriginPart = zip(range(len(OriginPart)), OriginPart)
        comp_mass = lambda x: len(x[1])
        OriginPart.sort(key=comp_mass)
        repo = set()
        NewPart = []
        for i, p in OriginPart:
            if len(p) > 1:
                p = set(p)
                p.difference_update(repo)
                repo.update(p)
            NewPart.append((i, list(p)))
        comp_id = lambda x: x[0]
        NewPart.sort(key=comp_id)
        NewPart = [x[1] for x in NewPart]
        return NewPart

    def LoadHostHalos(self, isnap=-1, selection=None):
        """Returns spatial properties of FoF groups for a snapshot.

        Arguments:
            isnap (int): snapshot number
            selection (tuple): (default = None) selection query being passed on
            to HBT file

        Returns:
            (numpy.ndarray) host haloes matching the selection
        """
        hosthalos = []
        offset = 0
        trans_index = False

        if selection is None:
            selection = np.s_[:]
        else:
            trans_index = isinstance(selection, numbers.Integral)

        if type(selection) is list:
            selection = tuple(selection)

        with self.OpenFile(isnap, filetype="HaloSize") as hostfile:
            nsub = hostfile["HostHalos"].shape[0]
            if trans_index:
                if offset + nsub > selection:
                    try:
                        hosthalos.append(
                            hostfile["HostHalos"][selection - offset])
                    except ValueError:
                        pass
                offset += nsub
            else:
                try:
                    hosthalos.append(hostfile["HostHalos"][selection])
                except ValueError:
                    pass

        if len(hosthalos):
            hosthalos = np.hstack(hosthalos)
        else:
            hosthalos = np.array(hosthalos)

        return hosthalos

    def GetHostHalo(self, HostHaloId, isnap=-1):
        """Returns spatial information of a specific FoF group.

        Arguments:
            HostHaloId (int): host halo id
            isnap (int): snapshot number

        Returns:
            (numpy.ndarray) single FoF group in a structured array
        """

        return self.LoadHostHalos(isnap, HostHaloId)[0]

    def GetSubsOfHost(self, HostHaloId, isnap=-1):
        """Loads all subhaloes belonging to a host halo.

        Uses information from ``Membership/GroupedTrackIds``.

        Arguments:
            HostHaloId (int): row number
            isnap (int): (default = -1)

        Returns:
            (numpy.ndarray): array of ``TrackIds``, or empty array if FoF does
                not host subhalos
        """
        with self.OpenFile(isnap) as subfile:
            trackIds = subfile['Membership/GroupedTrackIds'][HostHaloId]

        try:
            result = np.hstack(
                [self.GetSub(trackId, isnap) for trackId in trackIds])
        except ValueError:
            result = np.array([])

        return result

    def GetHostMergerTree(self, HostHaloId, isnap=-1, file=sys.stdout):
        """Builds a FoF merger tree starting at a host halo ID.

        Prints a graph which can be consumed by [Dot](graphviz.org] software.

        Arguments:
            HostHaloId (int): a halo at which to root a tree
            isnap (int): (default = -1) snapshot number
            file (File): (default=sys.stdout) default output for a Dot diagram

        Returns:
            (list): multiply embedded list with a tree, i.e. ``[1, [2, 3]]``
        """

        progenitors = self.GetProgenitorHaloes(HostHaloId, isnap)

        logging.debug("Halo %d at %d with %d progenitor(s)" %
                      (HostHaloId, isnap, progenitors.size))

        file.write("\t\"%d_%d\" [label=\"%d, %d\"]\n"%\
            (isnap, HostHaloId, isnap, HostHaloId))
        for progenitor in progenitors:
            file.write("\t\"%d_%d\" -- \"%d_%d\"\n" % (isnap, HostHaloId,
                                                       isnap - 1, progenitor))

        return [(HostHaloId, isnap), [] if progenitors.size == 0 else [
            self.GetHostMergerTree(progenitor, isnap - 1, file)
            for progenitor in progenitors
        ]]

    def GetProgenitorHaloes(self, HostHaloId, isnap=-1):
        """Returns progenitors of a given host halo.

        Given a HostHaloId and a snapshot number, it finds all subhaloes within
        the host, follows their track ids to the previous snapshot, checks if
        they are the main / central subhalo (by ``Rank`` values), finds those
        subhaloes' hosts and keeps the unique IDs.  Filtering by ``Rank`` is
        necessary to remove fly-by and re-merger events.

        Arguments:
            HostHaloId (int): host halo id
            isnap (int): (default=-1) snapshot number

        Returns:
            (numpy.ndarray): array of host halo ids of progenitors
        """

        try:
            subhaloes = self.GetSubsOfHost(HostHaloId, isnap=isnap)['TrackId']
            progenitors = np.unique([
                self.GetSub(subhalo, isnap - 1)[0]['HostHaloId']
                for subhalo in subhaloes
                if self.GetSub(subhalo, isnap - 1)[0]['Rank'] == 0
            ])
        except:
            progenitors = np.array([])

        return progenitors


if __name__ == '__main__':
    import timeit
    #apostle=HBTReader('../configs/Apostle_S1_LR.conf')
    apostle = HBTReader(
        '/cosma/home/jvbq85/data/HBT/data/apostle/S1_LR/subcat/VER1.8.1.param')
    #apostle=HBTReader('/cosma/home/jvbq85/data/HBT/data/MilliMill/subcat2_full/VER1.8.1.param')
    print(timeit.timeit(
        "[apostle.LoadSubhalos(i, 1) for i in range(10,apostle.MaxSnap)]",
        setup="from __main__ import apostle",
        number=1))
    #print(timeit.timeit("[apostle.LoadSubhalos(i, np.s_['Nbound','Rank']) for i in range(10,apostle.MaxSnap)]", setup="from __main__ import apostle,np", number=1))
    print(timeit.timeit(
        "[apostle.LoadSubhalos(i, 'Nbound') for i in range(10,apostle.MaxSnap)]",
        setup="from __main__ import apostle",
        number=1))
    print(timeit.timeit(
        "apostle.LoadSubhalos(-1, ('Nbound','Rank'))",
        setup="from __main__ import apostle",
        number=100))
    print(timeit.timeit(
        "[apostle.LoadSubhalos(i) for i in range(10,apostle.MaxSnap)]",
        setup="from __main__ import apostle",
        number=1))
    print(timeit.timeit(
        "apostle.GetTrack(12)", setup="from __main__ import apostle",
        number=1))
    print(timeit.timeit(
        "apostle.GetTrack(103)",
        setup="from __main__ import apostle",
        number=1))

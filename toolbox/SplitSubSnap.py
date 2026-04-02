# Copyright (C) 2007-2026 Jiaxin Han and contributors
# SPDX-License-Identifier: LGPL-3.0-or-later

''' Split each SubSnap file into SubTab and SubPart files
'''
import os,sys
import h5py

# rootdir=sys.argv[1]
run = sys.argv[1]
isnap = int(sys.argv[2])

rootdir = '/public/home/ackxyzh5n2/data/Jiutian/'+run
# rootdir=rootdir+'/subcat'

def SplitSnap(rootdir, isnap):
    print('snap', isnap)
    indir = rootdir+'/subcat-new/%03d/' % isnap
    outdir = rootdir+'/subcat-split/%03d/' % isnap
    os.makedirs(outdir, exist_ok=True)

    ifile = 0
    infilename = indir+'SubSnap_%03d.%d.hdf5' % (isnap, ifile)
    if not os.path.exists(infilename):
        print("Error:", infilename, "does not exist")
        return
    infile = h5py.File(infilename, 'r')
    nfiles = infile['NumberOfFiles'][0]
    infile.close()
    #print(nfiles)
    
    for ifile in range(nfiles):
        print(ifile, end=', ')
        infilename = indir+'SubSnap_%03d.%d.hdf5' % (isnap, ifile)
        infile = h5py.File(infilename, 'r')

        tabfilename = outdir+'SubTab_%03d.%d.hdf5' % (isnap, ifile)
        partfilename = outdir+'SubParticle_%03d.%d.hdf5' % (isnap, ifile)
        tabfile = h5py.File(tabfilename, 'w')
        partfile = h5py.File(partfilename, 'w')

        for item in infile.keys():
            if item == 'SubhaloParticles':  # |item=='ParticleProperties'
                infile.copy(item, partfile)
            else:
                infile.copy(item, tabfile)
        assert(tabfile['Subhalos'].shape == infile['Subhalos'].shape)
        assert(partfile['SubhaloParticles'].shape
               == infile['Subhalos'].shape)
        tabfile.close()
        partfile.close()
        infile.close()


if __name__ == '__main__':
    SplitSnap(rootdir, isnap)

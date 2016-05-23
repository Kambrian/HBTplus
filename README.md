new implementation of HBT in C++ . This is the OpenMP version for share memory machines. Check the [`MPI`](https://github.com/Kambrian/HBT2/tree/MPI) branch for a MPI version.

## Prerequisites

- a `c++` compiler with `c++11` support (e.g., `gcc 4.6.3` above)
- [HDF5](https://www.hdfgroup.org/) C library (1.8.0 and above)

### Optional dependence
- GNU Scientific Library [(GSL)](http://www.gnu.org/software/gsl/). 
Only needed if you want to output the shapes and orientations of subhaloes. To enable or disable GSL support, uncomment or comment out the GSL block in `Makefile.inc` (especially the `-DHAS_GSL` line). When enabled, HBT will do eigenvalue decomposition of the inertial tensor of each subhalo, and output the eigenvalue and eigenvectors describing the shape and direction of the subhalo. Without GSL, only the inertial tensors will be output.

## Compile
To produce single-precision `HBT` (internal datatypes are 4byte int and 4byte float), do

	make

. If you need double-precision, do

    make HBTdouble
    
## Run
 
    ./HBT configs/Example.conf [snapshotstart] [snapshotend]

Check `configs/Example.conf` for a sample parameter file.

If `snapshotend` is omitted, only process `snapshotstart`. If `snapshotstart` is also omitted, will run from `MinSnapshotIndex` (default=0) to `MaxSnapshotIndex` (specified in config file).

To submit to a batch queue, check `HBTjob.bsub`

## Output
The outputs are in HDF5 format, which can be viewed with [HDFView](https://www.hdfgroup.org/products/java/hdfview/index.html) or any other HDF tools. In python, you can use [h5py](https://pypi.python.org/pypi/h5py) to read them directly. There are two types of files in the output:
  
- SubSnap_*.hdf5: the subhalo catalogues.
- SrcSnap_*.hdf5: source subhalo catalogues, only used for restarting HBT from an intermediate snapshot. Can be safely removed once the run has finished.

Besides, the VER*.param records the version number of HBT used, as well as the parameter values used.

Each subhalo is labelled by a unique `TrackId`, which is fixed throughout its evolution history. So doing merger tree with HBT is straightforward: the progenitor/descendent of a subhalo at another snapshot is simply the subhalo labelled by the same `TrackId` at that time. The host halo of each subhalo is given by `HostId`, which is the index of the host halo in the order stored in the corresponding (FoF) halo catalogue. To facilitate fast retrieval of all the subhaloes in each host halo, the `/Membership/GroupedTrackIds` dataset in the file stores the list of subhaloes in each group (Note this is only available for the OpenMP version of HBT2). `Nbound` gives the number of bound particles in the subhalo. `Nbound=1` means the subhalo has been disrupted, so that only the most-bound particle is still tracked. The other properties should be self-explainatory.

Notes on Peebles and Bullock spin parameters: these parameters are vaguely defined due to the ambiguity/lack of standard in the mass, radius, and energy of a subhalo. Use them with caution. If possible, use the `SpecificAngularMomentum` instead of the spin parameters.

For the Hydrodynamical version of HBT, there could be objects with `Nbound=0` and an empty particle list. this means the track is lost due to all its particles consumed by a BH.

## Reference
For now, please cite the original [HBT paper](http://adsabs.harvard.edu/abs/2012MNRAS.427.2437H) if you use HBT in your work. We will soon have another paper coming out describing the new implementation here.

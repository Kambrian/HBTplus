new implementation of HBT in C++ . This is the hybrid MPI/OpenMP parallelized version. Check the [master](https://github.com/Kambrian/HBT2/tree/master) branch for a pure OpenMP version.

## Prerequisites

- a `c++` compiler with `c++11` support (e.g., `gcc 4.6.3` above)
- [HDF5](https://www.hdfgroup.org/) C library (1.8.0 and above)
- MPI library

### Optional dependence
- GNU Scientific Library [(GSL)](http://www.gnu.org/software/gsl/). 
Only needed if you want to output the shapes and orientations of subhaloes. To enable or disable GSL support, uncomment or comment out the GSL block in `Makefile.inc` (especially the `-DHAS_GSL` line). When enabled, HBT will do eigenvalue decomposition of the inertial tensor of each subhalo, and output the eigenvalue and eigenvectors describing the shape and direction of the subhalo. Without GSL, only the inertial tensors will be output.

## Compile
To produce single-precision `HBT` (internal datatypes are 4byte int and 4byte float), do

	make

. If you need double-precision, do

    make HBTdouble
    
### Special Compiler Flags
Below are a few macros to further customize the behaviour of HBT. These flags can be switched on or off in `Makefile.inc` and `Makefile`. Check the `CXXFLAGS` lines for these macros.
  
- `-DENABLE_EXPERIMENTAL_PROPERTIES`: output the peebles and bullock spin parameters. These parameters are vaguely defined due to the ambiguity/lack of standard in the mass, radius, and energy of a subhalo. Use them with caution. If possible, use the `SpecificAngularMomentum` instead of the spin parameters.

 
## Run

    mpirun -np 2 ./HBT configs/Example.conf [snapshotstart] [snapshotend]

will run it with 2 mpi processes. 

Check `configs/Example.conf` for a sample parameter file.

If `snapshotend` is omitted, only process `snapshotstart`. If `snapshotstart` is also omitted, will run from `MinSnapshotIndex` (default=0) to `MaxSnapshotIndex` (specified in config file).

To submit to a batch queue, check `HBTjob_mpi.bsub`

## Output
The outputs are in HDF5 format, which can be viewed with [HDFView](https://www.hdfgroup.org/products/java/hdfview/index.html) or any other HDF tools. In python, you can use [h5py](https://pypi.python.org/pypi/h5py) to read them directly. There are two types of files in the output:
  
- SubSnap_*.hdf5: the subhalo catalogues.
- SrcSnap_*.hdf5: source subhalo catalogues, only used for restarting HBT from an intermediate snapshot. Can be safely removed once the run has finished.

Besides, the VER*.param records the version number of HBT used, as well as the parameter values used.

Each subhalo is labelled by a unique `TrackId`, which is fixed throughout its evolution history. So doing merger tree with HBT is straightforward: the progenitor/descendent of a subhalo at another snapshot is simply the subhalo labelled by the same `TrackId` at that time. The host halo of each subhalo is given by `HostId`, which is the index of the host halo in the order stored in the corresponding (FoF) halo catalogue. To facilitate fast retrieval of all the subhaloes in each host halo, the `/Membership/GroupedTrackIds` dataset in the file stores the list of subhaloes in each group (Note this is only available for the OpenMP version of HBT2). `Nbound` gives the number of bound particles in the subhalo. `Rank` gives the order of subhaloes inside the group if sorted according to `Nbound`, with `Rank=0` indicating the most-massive subhalo inside each group (i.e., the main/central subhalo).

Once a subhalo is stripped to below `MinNumPartOfSub` specified in the parameter file, HBT continues to track its most bound particle. This single-particle descendents then have `Nbound=1`, and represent the "orphan" galaxy population in the framework of semi-analytical models. These orphans are also listed as subhaloes.

The other properties should be self-explainatory.

If you have difficulty reading the structure array of subhaloes or the variable length particle list, you can use `toolbox/convertSubSnap.py` to convert it to basic hdf5 files contains only vanilla arrays.

## Reference
For now, please cite the original [HBT paper](http://adsabs.harvard.edu/abs/2012MNRAS.427.2437H) if you use HBT in your work. We will soon have another paper coming out describing the new implementation here.


## Notes for users migrating from `HBT` to `HBT2`
HBT and HBT2 have different algorithmic details. They are not expected to give identical results. 

HBT no longer uses `ProSubID`. Instead, each subhalo is labelled by a unique `TrackId`, which is fixed throughout its evolution history. The progenitor/descendent of a subhalo at another snapshot is simply the subhalo labelled by the same `TrackId` at that time. 

sub_hierarchy is not available in HBT2 (but you would rarely need it.)

The host halo of each subhalo is given by `HostHaloId`, which is the index of the host halo in the order stored in the corresponding (FoF) halo catalogue.  With this you can sort or search to find all the members of each host.

HBT2 no longer have splintters. HBT2 does not store fake haloes either, i.e., for haloes that are not bound, you won't be able to find any subhalo hosted by it in HBT2.
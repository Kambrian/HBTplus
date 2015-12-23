new implementation of HBT in C++ . This is the MPI version.

## Prerequisites

- a `c++` compiler with `c++11` support (e.g., `gcc 4.6.3` above)
- [HDF5](https://www.hdfgroup.org/) C library (1.8.0 and above)
- MPI library

Core part done. More post-processing to be added.

## Compile
To produce single-precision `HBT` (internal datatypes are 4byte int and 4byte float), do

	make

. If you need double-precision, do

    make HBTdouble
 
## Run

    mpirun -np 2 HBT configs/Example.conf [snapshotstart] [snapshotend]

will run it with 2 mpi processes. Set `OMP_NUM_THREADS` to adjust the number of openmp threads on each node.

Check `configs/Example.conf` for a sample parameter file.

If `snapshotend` is omitted, only process `snapshotstart`. If `snapshotstart` is also omitted, will run from `MinSnapshotIndex` (default=0) to `MaxSnapshotIndex` (specified in config file).

To submit to a batch queue, check `HBTjob_mpi.bsub`

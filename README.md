new implementation of HBT in C++ . This is the OpenMP version for share memory machines. Check the `MPI` branch for a MPI version.

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
    
### Special Compiler Flags
Below are a few macros to further customize the behaviour of HBT. These flags can be switched on or off in `Makefile.inc` and `Makefile`. Check the `CXXFLAGS` lines for these macros.
  
- `-DENABLE_EXPERIMENTAL_PROPERTIES`: output the peebles and bullock spin parameters. These parameters are vaguely defined due to the ambiguity/lack of standard in the mass, radius, and energy of a subhalo. Use them with caution. If possible, use the `SpecificAngularMomentum` instead of the spin parameters.

- `-DALLOW_BINARY_SYSTEM`: give special treatment to binary systems-- those resulting from major mergers so that there is not a well-defined central subhalo. With this macro defined, HBT will not define a central subhalo for these systems, but treat all the subhaloes as satellite subhaloes. The mass ratio of the merger to define such binary systems is specified by the parameter `BinaryMassRatioLimit`. If you do not understand what I am talking about here, you probably do not need to care about it. This macro is enabled if you

    make HBTmajormerger
    
instead of `make HBT`.


 
## Run
 
    ./HBT configs/Example.conf [snapshotstart] [snapshotend]

Check `configs/Example.conf` for a sample parameter file.

If `snapshotend` is omitted, only process `snapshotstart`. If `snapshotstart` is also omitted, will run from `MinSnapshotIndex` (default=0) to `MaxSnapshotIndex` (specified in config file).

To submit to a batch queue, check `HBTjob.bsub`

## Reference
For now, please cite the original [HBT paper](http://adsabs.harvard.edu/abs/2012MNRAS.427.2437H) if you use HBT in your work. We will soon have another paper coming out describing the new implementation here.

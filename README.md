# HBT+ 
Hierarchical-Bound-Tracing (HBT) is a subhalo finder and merger tree builder that works by tracking the formation and evolution of subhalos over time. The tracking approach allows it to avoid common pitfalls in almost every other subhalo finder/merger tree builder. 

HBT+ is the new implementation in C++. It is completely rewritten for improved efficiency, usability, distributed memory support and with more complete physics.

This is the hybrid MPI/OpenMP parallelized version. Check the [Hydro](https://github.com/Kambrian/HBT2/tree/Hydro) branch for a pure OpenMP version.

Documentation available on the [wiki](https://github.com/Kambrian/HBT2/wiki).

## License and Citation

HBT+ is licensed under the [GNU Lesser General Public License v3.0 (LGPLv3)](./COPYING.LESSER). 

If you use this code for your research, please cite the following papers:
* Han et al. (2018): [MNRAS, 474, 604](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474..604H)
* Han et al. (2012): [MNRAS, 427, 2437](https://ui.adsabs.harvard.edu/abs/2012MNRAS.427.2437H)

## Contact
Jiaxin Han (SJTU): jiaxin.han _at_ sjtu.edu.cn

#!/bin/bash
#SBATCH --job-name=MiniJiutianG4-2 # 作业名
#SBATCH --no-requeue
#SBATCH --partition=kshcnormal01           # cpu 队列
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4      #threads; each node 32cpus
#SBATCH --output=logs/%x.%j
#SBATCH --error=logs/%x.%j
##SBATCH --mail-user=
##SBATCH --mail-type=ALL

ulimit -s unlimited
source addmodule
#module load compiler/devtoolset/7.3.1 mpi/hpcx/2.11.0/gcc-7.3.1 mathlib/gsl mathlib/hdf5/1.12.0/gnu731-cxx-fortran
#module load mpi/intelmpi compiler/intel mathlib/gsl mathlib/hdf5
#~ source /etc/profile.d/modules.csh

# by default, openmp is not enabled in the MPI version of HBT+. If you want to enable it, uncomment the openmp flags in Makefile.inc and rebuild.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# specify the first and last snapshot to process. useful for restarting from a specific snapshot. if not specified, then snapstart=0 and snapend=MaxSnapshotIndex.
#snapstart=65
#snapend=127

NP=$SLURM_NTASKS
echo $NP tasks each $OMP_NUM_THREADS threads

mpirun -np $NP ./HBTi8 configs/${SLURM_JOB_NAME}.conf $snapstart $snapend

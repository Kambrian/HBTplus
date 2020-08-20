#!/bin/bash
#SBATCH --job-name=PMO3000.hybrid # 作业名
#SBATCH --no-requeue
#SBATCH --partition=debug           # cpu 队列
#SBATCH --ntasks=64                    # mpi tasks
#SBATCH --cpus-per-task=30      #threads; each node 64cpus, so 2 tasks
#SBATCH --output=logs/%x.%j
#SBATCH --error=logs/%x.%j
#ls /public/home/jiaxinhan/PMO3000/simu/snapdir_*
#ls /public/home/liming/data/Simu_n3072_b500_pl13/snapdir_* -d

ulimit -s unlimited
module load mpi/intelmpi compiler/intel mathlib/gsl mathlib/hdf5
#~ source /etc/profile.d/modules.csh

# by default, openmp is not enabled in the MPI version of HBT+. If you want to enable it, uncomment the openmp flags in Makefile.inc and rebuild.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# specify the first and last snapshot to process. useful for restarting from a specific snapshot. if not specified, then snapstart=0 and snapend=MaxSnapshotIndex. 
snapstart=59
snapend=63

NP=$SLURM_NTASKS
echo $NP tasks each $OMP_NUM_THREADS threads

mpirun -np $NP ./HBTi8.hybrid configs/${SLURM_JOB_NAME}.conf $snapstart $snapend

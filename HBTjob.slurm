#!/bin/bash -l

#SBATCH --ntasks 28
#SBATCH -J AqA2
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive
#SBATCH -t 72:00:00
#SBATCH --mail-type=END                          # notifications for job done & fail
#SBATCH --mail-user=jiaxin.han@sjtu.edu.cn

module purge
#load the modules used to build your program.
module load intel_comp hdf5/1.10.3 gsl

export OMP_NUM_THREADS=$SLURM_NTASKS

# specify the first and last snapshot to process. useful for restarting from a specific snapshot. if not specified, then snapstart=0 and snapend=MaxSnapshotIndex. 
snapstart=98
snapend=162

# Run the program
./HBTi8 configs/${SLURM_JOB_NAME}.conf $snapstart $snapend


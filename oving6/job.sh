#!/bin/bash

#PBS -N ss_poisson

# Allocate three nodes with 12 processors from the default resources n*p=36
#PBS -lnodes=3:ppn=12:default

# Expect to run up to 5 minutes
#PBS -lwalltime=01:00:00

# Memory per process
#PBS -lpmem=2000MB

# Run on the freecycle account
#PBS -A freecycle

# Run in the optimist queue by default
#PBS -q optimist

# Join stdout and stderr output to one file
#PBS -j oe

#PBS -m abe
#PBS -M simen@hgranlund.com

# Change directory to dir with the job script
cd ${PBS_O_WORKDIR}

# Load needed modules
module load intelcomp
module load openmpi/1.4.3-intel

# Set thread affinity
KMP_AFFINITY="granularity=fine,compact"


# Run with 8 MPI processes, each with 3 threads
OMP_NUM_THREADS=1 mpirun -npernode 12 oving6-poison-mpi 5 14
OMP_NUM_THREADS=2 mpirun -npernode 6 oving6-poison-mpi 5 14
OMP_NUM_THREADS=3 mpirun -npernode 4 oving6-poison-mpi 5 14
OMP_NUM_THREADS=4 mpirun -npernode 3 oving6-poison-mpi 5 14
OMP_NUM_THREADS=6 mpirun -npernode 2 oving6-poison-mpi 5 14
OMP_NUM_THREADS=12 mpirun -npernode 1 oving6-poison-mpi 5 14

#!/bin/bash

#SBATCH --time=0:4:0   # walltime
#SBATCH --nodes=8   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=16384M   # memory per CPU core
#SBATCH -J "Mandelbrot"   # job name
#SBATCH --mail-user=gseccles@gmail.com   # email address
#SBATCH --verbose

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
mpirun parallelMand


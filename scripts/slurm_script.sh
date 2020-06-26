#!/bin/bash
#
# job name:
#SBATCH -J fdtd
#
# number of jobs rounded % 8:
#SBATCH -n 64
#
mpirun -np 64 ./mpisolve

#!/bin/bash
#SBATCH -n 2 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-00:05 # Maximum execution time (D-HH:MM)
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=512 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname

mpirun -np 2 ./ip-catchall mhInput.inp

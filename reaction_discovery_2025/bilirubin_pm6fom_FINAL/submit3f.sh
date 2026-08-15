#!/usr/bin/bash
#SBATCH -J min           # Job name
#SBATCH -o sout      # Name of stdout output file
#SBATCH -e serror       # Name of stderr error file
#SBATCH -p gpu-a100            # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
##SBATCH --gpus-per-node=3  # Number of GPU(s) per node
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --cpus-per-task 4
#SBATCH -A CHE23043       # Project/Allocation name (req'd if you have more than 1)
##SBATCH --mail-user=ankit.pandey@ttu.edu

cd $SLURM_SUBMIT_DIR

cd secondary6
./submit.sh &
cd ..
sleep 1

cd secondary16
./submit.sh &
cd ..
sleep 1

cd secondary26
./submit.sh &
cd ..

wait


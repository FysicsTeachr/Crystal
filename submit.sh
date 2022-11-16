#!/bin/bash
#SBATCH -J ch2nh
#SBATCH -o ch2nh.o
#SBATCH -e ch2nh.e
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p small
#SBATCH -t 01:00:00

export OMP_NUM_THREADS=56

module purge
module use /work2/01114/jfonner/frontera/modulefiles
module load gaussian/16rC.01

./a.out &
sleep 10
for i in {1..100}
  do
    g16 gin.gjf
    sleep 10
  done



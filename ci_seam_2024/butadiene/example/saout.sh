#!/usr/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=matador 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH -J terachem
#SBATCH -o terachem.o%j
#SBATCH -e terachem.e%j

rm shell*
rm set0*
rm foundp*
./a.out


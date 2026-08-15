#!/usr/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=matador 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH -J ts2
#SBATCH -o terachem.o%j
#SBATCH -e terachem.e%j

cd $SLURM_SUBMIT_DIR

source ~/geodesic-env/bin/activate
python -m geodesic_interpolate --nimages 15 --output interpolated.xyz 12.xyz

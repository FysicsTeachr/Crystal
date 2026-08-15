#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name GLU
#SBATCH --output=test2/slurm-out
#SBATCH --error=test2/slurm-error
#SBATCH --mem=128GB

cd $SLURM_SUBMIT_DIR

 source ~/geodesic-env/bin/activate
geodesic_interpolate 29-51_dft.xyz --output interpolated_29-51b.xyz --nimages 15 


#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name GLU
#SBATCH --output=test2/slurm-out
#SBATCH --error=test2/slurm-error
#SBATCH --mem=4GB

#cd $SLURM_SUBMIT_DIR

#chemsh pun2xyz.chm
for file in ../yesfolder/*.pun; do
    # Extract the base name without extension
    base_name=$(basename "$file" .pun)

    # Create the command to convert .pun to .xyz
    echo "write_xyz file=${base_name}.xyz coords=${file}" > pun2xyz.chm

    # Run chemsh with the generated .chm file
    chemsh pun2xyz.chm
done

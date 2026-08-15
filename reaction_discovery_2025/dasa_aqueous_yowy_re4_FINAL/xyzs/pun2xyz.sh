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
# Initialize an empty list
commands=""

# Loop through all .pun files in the folder
for file in ../yesfolder/*.pun; do
    # Extract the base name without extension
    base_name=$(basename "$file" .pun)

    # Add the command to the list
    commands+="write_xyz file=${base_name}.xyz coords=${file}\n"
done

# Write all commands to the pun2xyz.chm file in one go
echo -e "$commands" > pun2xyz.chm

# Run ChemShell once with the generated .chm file
chemsh pun2xyz.chm


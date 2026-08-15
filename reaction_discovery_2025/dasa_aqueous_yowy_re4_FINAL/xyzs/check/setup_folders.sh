#!/bin/bash
# Ensure the folders directory exists
mkdir -p folders
# Iterate through all .xyz files in ../yesfolder
for file in ../*.xyz; do
    # Extract the base filename without extension
    filename=$(basename -- "$file")
    folder_name="${filename%.*}"
    # Create a new folder inside the folders directory
    mkdir -p "folders/$folder_name"
    # Copy the .xyz file to the corresponding folder
    cp "$file" "folders/$folder_name/$folder_name.xyz"
    # Copy the a.out file into the folder
    cp a.out "folders/$folder_name"
    # Create the submit.sh script inside the folder
    cat > "folders/$folder_name/submit.sh" <<EOL
#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name minimize
#SBATCH --output=test/slurm-out
#SBATCH --error=test/slurm-error
#SBATCH --mem=4GB
cd \$SLURM_SUBMIT_DIR
source ~/.bashrc
source /home/pan60047/anaconda3/etc/profile.d/conda.sh
conda activate python3p8

cp $folder_name.xyz pregeom.xyz
./a.out
geometric-optimize --engine mopac --mopacexe /home/pan60047/mopacpi-master/mopacpi.x gradient.dat --maxiter 1000
EOL
done

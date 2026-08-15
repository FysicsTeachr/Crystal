#!/usr/bin/bash
#SBATCH -J min           # Job name
#SBATCH -o sout      # Name of stdout output file
#SBATCH -e serror       # Name of stderr error file
#SBATCH -p gpu-a100-small            # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
##SBATCH --gpus-per-node=1         # Number of GPU(s) per node
#SBATCH -n 1              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 2:00:00        # Run time (hh:mm:ss)
#SBATCH --cpus-per-task 4
#SBATCH -A CHE23021       # Project/Allocation name (req'd if you have more than 1)
##SBATCH --mail-user=ankit.pandey@ttu.edu

cd $SLURM_SUBMIT_DIR

conda init
conda activate py39
export TeraChem=/work/09717/ruibinliang/ls6/terachem_06_2023/terachem/TeraChem
export NBOEXE=/work/09717/ruibinliang/ls6/terachem_06_2023/terachem/TeraChem/bin/nbo6.i4.exe
export LD_LIBRARY_PATH=/work/09717/ruibinliang/ls6/terachem_06_2023/terachem/TeraChem/lib:$LD_LIBRARY_PATH
export PATH=/work/09717/ruibinliang/ls6/terachem_06_2023/terachem/TeraChem/bin:$PATH
export OpenMM=/work/04045/ankipand/ls6/anaconda3/
export OPENMM_PLUGIN_DIR=$OpenMM/lib/plugins

/work/09717/ruibinliang/ls6/terachem_06_2023/terachem/TeraChem/bin/terachem move_1.in  > ttc.out


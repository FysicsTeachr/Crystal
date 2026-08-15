#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name GLU 
#SBATCH --output=slurm-o
#SBATCH --error=slurm-e
#SBATCH --mem=4GB

cd $SLURM_SUBMIT_DIR
source ~/.bashrc

# Turbomole
source ~/source_turbomole.sh 
export PARA_ARCH=SMP
export PARNODES=4
export PATH=$TURBODIR/bin/x86_64-unknown-linux-gnu_smp/:$PATH

# ORCA
export ORCADIR=/lustre/work/rliang/orca_5_0_3_linux_x86-64_shared_openmpi411/
export PATH=$ORCADIR:$PATH
export LD_LIBRARY_PATH=$ORCADIR:$LD_LIBRARY_PATH

# ChemShell
export PATH=/lustre/work/rliang/chemsh-3.7.1/scripts:$PATH/

# xTB
export OMP_STACKSIZE=4G
ulimit -s unlimited
source /home/gujulian/PROGRAMS/xtb-dist/share/xtb/config_env.bash

for file in ../yesfolder/*.pun; do
        filename=$(basename -- "$file" .pun)
	cp "$file" example.pun
	chemsh example.chm
        mv optimized_1.pun outputs/optimized$filename.pun
	mv slurm-o outputs/o$filename
done


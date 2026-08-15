#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name pm7s
#SBATCH --output=test/slurm-out
#SBATCH --error=test/slurm-error
#SBATCH --mem=4GB
cd $SLURM_SUBMIT_DIR
for i in {1288..6448}; do
	cd $i
	sbatch submit.sh
	cd ..
done

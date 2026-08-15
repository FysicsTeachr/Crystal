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
source ~/.bashrc
source /home/pan60047/anaconda3/etc/profile.d/conda.sh
conda activate python3p8

##mv *.xyz pregeom.xyz
cp ../convert.out .
./convert.out
geometric-optimize --engine mopac --mopacexe /home/pan60047/mopacpi-master/mopacpi.x gradient.dat --converge energy 1.0e-5 grms 6.0e-4 gmax 9.0e-4 drms 2.4e-3 dmax 3.6e-3 --maxiter 9999


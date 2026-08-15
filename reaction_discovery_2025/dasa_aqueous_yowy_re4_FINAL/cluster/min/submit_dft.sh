#!/usr/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=matador 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH -J omega_min
#SBATCH -o terachem.o%j
#SBATCH -e terachem.e%j

cd $SLURM_SUBMIT_DIR
source ~/.bashrc
conda activate python3p8

source /home/rliang/intel/parallel_studio_xe_2019.5.075/bin/psxevars.sh
export SPACK_ROOT=/home/pan60047/spack
export PATH=$SPACK_ROOT/bin:$PATH
. /home/pan60047/spack/share/spack/setup-env.sh
export TeraChem=/lustre/work/rliang/terachem-tip/production/build
export LD_LIBRARY_PATH=$TeraChem/lib:/home/pan60047/anaconda3/envs/python3/lib:$LD_LIBRARY_PATH
module load protobuf-3.12.2-gcc-9.3.0-x4lhuqp
module load gcc/9.3.0
module load cuda-11.0.2-gcc-9.3.0-aji4ewc
module load libmatheval-1.1.11-gcc-9.3.0-3xuodfl
export OpenMM=/home/pan60047/anaconda3/envs/python3p8
export OPENMM_PLUGIN_DIR=$OpenMM/lib/plugins
export LD_LIBRARY_PATH=/usr/local/gcc/9.3.0/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/pan60047/anaconda3/envs/python3p8/lib:$LD_LIBRARY_PATH

$TeraChem/bin/terachem min_dft.in > dft.out 

#!/usr/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=matador 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH -J terachem
#SBATCH -o terachem.o%j
#SBATCH -e terachem.e%j

cd $SLURM_SUBMIT_DIR
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

rm foundli*
rm tgeom.xyz
rm ttc.out
rm ttc2.out

while true; do
    # Read the first line of the file
    first_line=$(head -n 1 "../start.txt")

    if [ "$first_line" == "1" ]; then
        # The first line is '1', break the loop
        echo "First line is '1', proceeding with the script..."
        break
    else
        # The first line is not '1', sleep for 2 seconds
        #echo "First line is not '1', sleeping for 2 seconds..."
        sleep 2
    fi
done

sleep 10

for i in {1..999999}
  do
    rm -r scr*
    rm terachem.*
    ./a.out
    cancl=$(head -n 1 cancel.txt)
    if [[ $cancl == *n* ]]; then
      $TeraChem/bin/terachem ttc.in > ttc.out
    fi
    ./b.out
  done


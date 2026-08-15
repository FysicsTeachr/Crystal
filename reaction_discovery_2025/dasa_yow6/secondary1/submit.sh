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

rm foundli*
rm tgeom.xyz
rm ttc.out
rm ttc2.out
rm terachem.*

terminate_after_timeout() {
  sleep 180 # 10 minutes in seconds
  if ps -p $1 > /dev/null; then
    kill -9 $1
    echo "y" > cancel.txt
  fi
}

sleep 10

for i in {1..999999}
do
  rm -r scr*
  rm terachem.*

  ./a.out
  python python.py
  python dissociation.py  
  if grep -q "yes" dissociated.txt; then
    read -r A B < coord.txt
    echo "Job terminated" > ttc.out
    mkdir -p ../dissociated
    cp tgeom.xyz ../dissociated/tgeom${A}.xyz
  else
    # Run TeraChem with monitoring for timeout
    $TeraChem/bin/terachem ttc.in > ttc.out &
    pid=$!
    terminate_after_timeout $pid &
    wait $pid
    exit_status=$?
    if [[ $exit_status -ne 0 ]]; then
      read -r A B < coord.txt
      cp ttc.out ../failed/ttc${A}.out
      cp tgeom.xyz ../failed/tgeom${A}.xyz
    fi
  fi
  
  ./b.out

  cancl=$(head -n 1 cancel.txt)
  if [[ $cancl == *n* ]]; then
    read -r A B < coord.txt
    echo ${A} > temp.txt
    tail -n 44 scr.tgeom/optim.xyz >> temp.txt

    # Call Python script to handle writing to allfiles.txt
    python temp.py temp.txt

    echo ${A} >> ../filesinfo2.txt # Append the marker to filesinfo2.txt
  fi
  echo "d" > status.txt
done


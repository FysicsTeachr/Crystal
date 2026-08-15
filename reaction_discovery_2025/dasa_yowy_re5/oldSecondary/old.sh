#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name re_3
#SBATCH --output=test2/slurm-out
#SBATCH --error=test2/slurm-error
#SBATCH --mem=4GB

cd $SLURM_SUBMIT_DIR
##source ~/.bashrc
source /home/pan60047/anaconda3/etc/profile.d/conda.sh
conda activate python3p8

# ORCA
export ORCADIR=/lustre/work/rliang/orca_5_0_3_linux_x86-64_shared_openmpi411/
export PATH=$ORCADIR:$PATH
export LD_LIBRARY_PATH=$ORCADIR:$LD_LIBRARY_PATH

# ChemShell
export PATH=/lustre/work/rliang/chemsh-3.7.1/scripts:$PATH/

# xTB
ulimit -s unlimited
source /home/gujulian/PROGRAMS/xtb-dist/share/xtb/config_env.bash

# Function to terminate a process after a timeout
terminate_after_timeout() {
  sleep 36000  # Timeout of 5 minutes (300 seconds)
  if ps -p $1 > /dev/null; then
    kill -9 $1
    echo "Timeout occurred, moving 1.pun and replacing slurm-out with sm.out"
    cp ../../1.pun ../../allfolder/${A}.pun
    cp sm.out slurm-out
  fi
}

ml intel
#ifort s1.f90
#ifort s2.f90 -o b.out
rm foundling*
sleep 10

for i in {1..999999}
do
  ./a.out
  read -r A B < coord.txt
  python makechm.py
  cd test2
  rm optimized.pun
  cp ../pregeom.pun .
  cp ../move.chm .
  > slurm-out
  
  # Run chemsh with timeout monitoring
  chemsh move.chm &
  pid=$!
  terminate_after_timeout $pid &  # Start timeout function in background
  wait $pid                      # Wait for chemsh process to finish
  exit_status=$?

  if [[ $exit_status -eq 0 ]]; then
    cp optimized_3.pun optimized.pun
    mv optimized.pun ../../allfolder/${A}.pun
  fi

    cp slurm-out allslurms/${A}x${B}-out
  
  cd ..
  ./b.out
  echo ${A} >> ../filesinfo2.txt  # Append the marker to filesinfo2.txt
  echo "d" > status.txt
done


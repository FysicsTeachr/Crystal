#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name re_4
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

terminate_if_stuck() {
  pid=$1
  slurm_file="slurm-out"
  timeout=1200  # Timeout in seconds
  check_interval=30  # Interval to check the file update (in seconds)
  max_idle_time=90  # Maximum idle time for output updates (in seconds)

  start_time=$(date +%s)
  last_update_time=$(date +%s)

  while true; do
    sleep $check_interval

    # Check if process is still running
    if ! ps -p $pid > /dev/null; then
      return  # Process has finished
    fi

    # Check last modification time of slurm-out
    if [ -e "$slurm_file" ]; then
      current_update_time=$(stat -c %Y "$slurm_file")
      if [ $current_update_time -gt $last_update_time ]; then
        last_update_time=$current_update_time  # Update last update time
      fi

      # Check for energy condition
      if grep -q "Final converged energy:" "$slurm_file"; then
        energy=$(grep "Final converged energy:" "$slurm_file" | tail -n 1 | awk '{print $4}')
        # Check if energy is in the specified range
        if (( $(echo "$energy > -71.2" | bc -l) && $(echo "$energy < -70" | bc -l) )); then
          echo "Energy condition met: $energy. Terminating process."
          kill -9 $pid
          echo "Moving 1.pun and replacing slurm-out with sm.out"
          cp ../../1.pun ../../allfolder/${A}.pun
          cp sm.out slurm-out
          return
        fi
      fi
    fi

    # Check for timeout or idle condition
    current_time=$(date +%s)
    elapsed_time=$((current_time - start_time))
    idle_time=$((current_time - last_update_time))

    if [ $elapsed_time -ge $timeout ] || [ $idle_time -ge $max_idle_time ]; then
      kill -9 $pid
      echo "Timeout or idle detected, terminating process."
      echo "Moving 1.pun and replacing slurm-out with sm.out"
      cp ../../1.pun ../../allfolder/${A}.pun
      cp sm.out slurm-out
      return
    fi
  done
}


for i in {1..999999}
do
  python s1.py  # Execute s1.py as in the original script
  read -r A B < coord.txt
  python makechm.py
  cd test2
  rm optimized.pun
  cp ../pregeom.pun .
  cp ../move.chm .
  > slurm-out

  # Run chemsh with monitoring
  chemsh move.chm &
  pid=$!
  terminate_if_stuck $pid &  # Start monitoring function in the background
  wait $pid                 # Wait for chemsh process to finish
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


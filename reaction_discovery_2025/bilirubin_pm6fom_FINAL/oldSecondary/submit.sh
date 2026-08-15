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


terminate_after_timeout() {
    local pid=$1
    local timeout=3600
    ( sleep $timeout && kill -9 $pid > /dev/null 2>&1 && echo "Process $pid killed after timeout." ) &
}

rm foundling*
rm values.txt
sleep 2

num_steps=3     ##     <-------------------------- steps

for i in {1..999999}
do
#  ./a.out
  python s1.py
  read -r A B < coord.txt
  python maketera.py
  cd test
  cp ../pregeom.xyz .
  ./convert.out

for step in $(seq 1 $num_steps)
do
   cp ../constraints_${step}.txt constraints.txt
   geometric-optimize --engine mopac --mopacexe /home/pan60047/mopacpi-master/mopacpi.x gradient.dat constraints.txt --maxiter 60 --converge energy 2.0e-5 grms 6.0e-3 gmax 9.0e-3 drms 2.4e-2 dmax 3.6e-2 &
   pid=$!
   terminate_after_timeout $pid &
   wait $pid
   exit_status=$?

   echo "geometric-optimize exit status: $exit_status"
   if [ $exit_status -ne 0 ]; then
       echo "geometric-optimize failed with status $exit_status"
   fi

   echo "Executing ./c.out"
   ./c.out
   exit_status_c=$?
   echo "./c.out exit status: $exit_status_c"
   if [ $exit_status_c -ne 0 ]; then
       echo "./c.out failed with status $exit_status_c"
   fi

   if [ ! -s proceed.txt ]; then
       echo "proceed.txt is missing or empty!"
   fi

   read -r Cs < proceed.txt
   echo "Value of Cs: $Cs"
   if [[ $Cs == "n" ]]; then
       echo "Cs is 'n', breaking loop."
       mv gradient.log allslurms/${A}x${B}x${step}x.out
       break
   elif [[ $Cs == "y" ]]; then
       echo "Cs is 'y', processing pregeom.xyz."
       tail -n 81 gradient_optim.xyz > pregeom.xyz
       ./convert.out
       echo " converted from output to next input"

   else
       echo "Unexpected value in proceed.txt: $Cs"
   fi

   mv gradient.log allslurms/${A}x${B}x${step}.out
   cp gradient_optim.xyz allslurms/${A}x${B}x${step}.xyz

done

#  rm -r scr*
  mv pregeom.xyz ../../allfolder/${A}.xyz
  if [ $? -ne 0 ]; then
    cp ../../1.xyz ../../allfolder/${A}.xyz
  fi
##  cp ttc.out allslurms/${A}x${B}-out
  cd ..
  ./b.out
  rm values.txt


  echo ${A} >> ../filesinfo2.txt  
  echo "d" > status.txt
done


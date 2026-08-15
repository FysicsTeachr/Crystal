#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name GLU
#SBATCH --output=test2/slurm-out
#SBATCH --error=test2/slurm-error
#SBATCH --mem=4GB
cd $SLURM_SUBMIT_DIR
source ~/.bashrc

terminate_after_timeout() {
  sleep 9999
  if ps -p $1 > /dev/null; then
    kill -9 $1
  fi
}

rm foundling*
rm values.txt
sleep 2

num_steps=3     ##     <-------------------------- steps

for i in {1..999999}
do
  ./a.out
  read -r A B < coord.txt
  python maketera.py
  cd test
  cp ../pregeom.xyz .
  ./convert.out

  for step in $(seq 1 $num_steps)
  do
     cp ../constraints_${step}.txt constraints.txt
     geometric-optimize --engine mopac --mopacexe /home/pan60047/mopacpi-master/mopacpi.x gradient.dat constraints.txt &
     pid=$!
     terminate_after_timeout $pid &  
     wait $pid  
     exit_status=$?
       ./c.out
       read -r Cs < proceed.txt
     if [[ $Cs == "n" ]]; then 
       break
     fi
     if [[ $Cs == "y" ]]; then
       tail -n 81 gradient_optim.xyz > pregeom.xyz          ##  <------------------------ previous scr filename
       ./convert.out
     fi
  done

  rm -r scr*
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


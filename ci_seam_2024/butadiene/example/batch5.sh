#!/usr/bin/bash

##./loadFoundlings.sh

cd secondary1
sbatch submit.sh
cd ..
cd secondary2
sbatch submit.sh
cd ..
cd secondary3
sbatch submit.sh
cd ..
cd secondary4
sbatch submit.sh
cd ..
cd secondary5
sbatch submit.sh
cd ..

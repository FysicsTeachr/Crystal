#!/usr/bin/bash

for i in {1..399}
do
    cd secondary$i
    sbatch submit.sh
    cd ..
done


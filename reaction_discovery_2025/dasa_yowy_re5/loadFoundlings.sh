#!/usr/bin/bash

# Define the variable xt
xt=399

# Remove old secondary directories
rm -r secondary*

# Loop from 1 to the value of xt
for i in $(seq 1 $xt)
do
    cp -R oldSecondary secondary$i

    # Modify the condition in s1.f90 using sed
#    sed -i "s/if(irr\/=0\.or\.current==0)/if(irr\/=0 .or. (current<$xt .and. current\/=$i))/g" secondary$i/s1.f90
done


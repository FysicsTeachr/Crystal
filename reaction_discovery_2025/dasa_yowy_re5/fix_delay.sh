#!/bin/bash

ml intel
# Loop through each value of i
for i in {1..399}; do  # Adjust the range {1..10} as needed
    file="secondary$i/s1.f90"

    # Check if the file exists before making changes
#    if [[ -f $file ]]; then
        # Use sed to find and replace the line
        sed -i "s/if(irr\/=0.or.current==0)then/if(irr\/=0.or.current==0 .or. current<$i)then/" "$file"
        echo "Updated $file"
#    else
#        echo "File $file not found"
#    fi
cd secondary$i
ifort s1.f90
cd ..

done

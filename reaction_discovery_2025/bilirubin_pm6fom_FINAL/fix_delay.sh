#!/bin/bash

# Loop through each value of i
for i in {1..399}; do  # Adjust the range {1..399} as needed
    file="secondary$i/s1.py"

    # Check if the file exists before making changes
    if [[ -f $file ]]; then
        # Use sed to replace the condition directly
        sed -i "s/if current == 0:/if current == 0 or current < $i:/" "$file"
        echo "Updated $file"
    else
        echo "File $file not found"
    fi

    # Check if the Python script runs without errors
#    cd secondary$i
 #   python3 -m py_compile s1.py
  #  cd ..
done


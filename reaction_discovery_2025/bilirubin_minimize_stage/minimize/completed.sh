#!/bin/bash

# Ensure the completed directory exists
mkdir -p completed

# Iterate through each subdirectory in folders/
for dir in folders/*/; do
    # Check if gradient.log exists in the current subdirectory
    if [ -f "$dir/gradient.log" ]; then
        # Get the last line of gradient.log
        last_line=$(tail -n 1 "$dir/gradient.log")
        
        # Check if the last line starts with "Time elapsed"
        if [[ $last_line == \ \ \ \ Time\ elapsed* ]]; then
            # Check if gradient_optim.xyz exists in the current subdirectory
            if [ -f "$dir/gradient_optim.xyz" ]; then
                # Extract the last 81 lines and save to completed/<subdir_name>.xyz
                subdir_name=$(basename "$dir")
                tail -n 81 "$dir/gradient_optim.xyz" > "completed/${subdir_name}.xyz"
                echo "Processed: completed/${subdir_name}.xyz"
            else
                echo "gradient_optim.xyz not found in $dir"
            fi
        fi
    else
        echo "gradient.log not found in $dir"
    fi
done


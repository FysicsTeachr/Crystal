#!/bin/bash

# Iterate through all directories inside folders
for dir in folders/*/; do
    # Change into the directory
    cd "$dir"
    
    # Submit the submit.sh script
    sbatch submit.sh
    
    # Return to the original directory
    cd - > /dev/null
done


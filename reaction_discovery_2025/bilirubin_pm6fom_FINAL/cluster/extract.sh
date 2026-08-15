#!/bin/bash

# Create or clear the output directory
OUTPUT_DIR="minimized_2"
if [ -d "$OUTPUT_DIR" ]; then
    rm -rf "$OUTPUT_DIR"/*
else
    mkdir "$OUTPUT_DIR"
fi

# Loop through folders from 1 to 6448
for i in {1..6448}; do
    # Check if the folder and gradient.log file exist
    FOLDER="../min/$i"
    LOG_FILE="$FOLDER/gradient.log"
    XYZ_FILE="$FOLDER/gradient_optim.xyz"
    
    if [ -f "$LOG_FILE" ] && tail -n 20 "$LOG_FILE" | grep -q "Converged! =D"; then
        # Extract XYZ data if the calculation was completed
        if [ -f "$XYZ_FILE" ]; then
            tail -n 81 "$XYZ_FILE" > "$OUTPUT_DIR/$i.xyz"
        else
            echo "Warning: XYZ file missing for folder $i."
        fi
    fi
done

echo "Extraction completed. Check the '$OUTPUT_DIR' directory for results."


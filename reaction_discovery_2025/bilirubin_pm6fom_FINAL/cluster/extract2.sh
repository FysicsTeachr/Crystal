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
    # Define folder paths
    FOLDER="../min/$i"
    LOG_FILE="$FOLDER/gradient.log"
    BACKUP_LOG_FILE="$FOLDER/gradient_0.log"
    XYZ_FILE="$FOLDER/gradient_optim.xyz"
    BACKUP_XYZ_FILE="$FOLDER/gradient0_optim.xyz"

    # Check if the primary log file exists and contains "Converged! =D"
    if [ -f "$LOG_FILE" ] && tail -n 20 "$LOG_FILE" | grep -q "Converged! =D"; then
        # Extract XYZ data from the primary XYZ file
        if [ -f "$XYZ_FILE" ]; then
            tail -n 81 "$XYZ_FILE" > "$OUTPUT_DIR/$i.xyz"
        else
            echo "Warning: XYZ file missing for folder $i."
        fi
    # If primary log doesn't contain the phrase, check the backup log
    elif [ -f "$BACKUP_LOG_FILE" ] && tail -n 20 "$BACKUP_LOG_FILE" | grep -q "Converged! =D"; then
        # Extract XYZ data from the backup XYZ file
        if [ -f "$BACKUP_XYZ_FILE" ]; then
            tail -n 81 "$BACKUP_XYZ_FILE" > "$OUTPUT_DIR/$i.xyz"
        else
            echo "Warning: Backup XYZ file missing for folder $i."
        fi
    fi
done

echo "Extraction completed. Check the '$OUTPUT_DIR' directory for results."


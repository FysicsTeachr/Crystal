#!/bin/bash

# Clear or create the output file
output_file="all.xyz"
> "$output_file"

# Find, sort, and batch concatenate
ls *.xyz 2>/dev/null | sort -t. -k1,1n | xargs cat >> "$output_file"

echo "Files concatenated into $output_file in numerical order."


#!/usr/bin/bash

input_file="unique_geometries.xyz"

# Extract IDs from the file into an array
ids=($(awk '/^Energy/ {print $NF}' "$input_file"))

for i in "${ids[@]}"; do
    mkdir -p "$i"
    cp "../../minimize/$i/optimized.pun" "$i/"
done


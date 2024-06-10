#!/bin/bash

# Define the prefix of your files
prefix="ref_amr2_nt200_tfin20"

# Initialize an empty array to hold the numeric parts of the filenames
declare -a indices

# Populate the array with the numeric parts extracted from filenames
for file in ${prefix}.*.hdf5; do
    index=$(echo "$file" | sed "s/${prefix}\.\(.*\)\.hdf5/\1/")
    indices+=("$index")
done

# Find the maximum index in the array
max_index=$(printf "%s\n" "${indices[@]}" | sort -n | tail -1)

# Loop through all the hdf5 files with the specified prefix again
for file in ${prefix}.*.hdf5; do
    # Extract the numeric part of the current file
    index=$(echo "$file" | sed "s/${prefix}\.\(.*\)\.hdf5/\1/")
    # Check if the file is not the first or the last one
    if [[ "$index" != "0000000" && "$index" != "$max_index" ]]; then
        # Remove the file
        rm "$file"
    fi
done

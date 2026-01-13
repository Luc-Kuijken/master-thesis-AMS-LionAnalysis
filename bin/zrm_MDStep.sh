#!/bin/bash

# Developed by L.C. Kuijken
# LAST UPDATED 26-07-25

search_path="${1-}"

# Check if directory exists
if [[ ! -d "$search_path" ]]; then
    echo "Directory does not exist: $search_path"
    exit 1
fi

# Find and delete files starting with 'MDStep'
find "$search_path" -type f -name 'MDStep*' -exec rm -v {} \;

echo "Finished deleting all MDStep* files in $search_path"